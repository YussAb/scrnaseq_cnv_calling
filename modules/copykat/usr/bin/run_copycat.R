#!/usr/bin/env Rscript

# Lightweight CopyKAT runner used by the Nextflow module. It translates the
# Nextflow process inputs into the argument names expected by copykat().

# Parse command-line arguments passed as --key value pairs. Arguments without a
# following value are treated as boolean flags.
parse_args <- function(args) {
    parsed <- list()
    i <- 1

    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--")) {
            stop("Unexpected positional argument: ", key, call. = FALSE)
        }

        name <- sub("^--", "", key)
        next_i <- i + 1
        if (next_i > length(args) || startsWith(args[[next_i]], "--")) {
            parsed[[name]] <- TRUE
            i <- i + 1
        } else {
            parsed[[name]] <- args[[next_i]]
            i <- i + 2
        }
    }

    parsed
}

# Fetch a required argument and stop with a clear message if it is missing.
required_arg <- function(args, name) {
    value <- args[[name]]
    if (is.null(value) || identical(value, "")) {
        stop("Missing required argument --", name, call. = FALSE)
    }
    value
}

# Small default-value helper: lhs %||% rhs returns rhs only when lhs is NULL.
`%||%` <- function(lhs, rhs) {
    if (is.null(lhs)) rhs else lhs
}

# Convert Nextflow-style boolean values to R logical values.
as_bool <- function(value) {
    tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

# Convert comma-separated parameters, such as ref_group_names, into a character
# vector. Empty values become character(0) so downstream checks are simple.
as_csv_vector <- function(value) {
    if (is.null(value) || isTRUE(value) || identical(value, "")) {
        return(character())
    }
    trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
}

# Keep boolean options in the TRUE/FALSE string form passed to copykat().
copykat_bool <- function(value) {
    if (as_bool(value)) "TRUE" else "FALSE"
}

# Load CopyKAT inside the process and fail early if the execution environment is
# missing the package.
load_copykat <- function() {
    if (!suppressPackageStartupMessages(require("copykat", character.only = TRUE))) {
        stop(
            "The R package 'copykat' is not installed. ",
            "Run with a container/conda environment that provides copykat.",
            call. = FALSE
        )
    }

}

# CopyKAT expects a numeric matrix with genes as row names and cells as columns.
# The input file is therefore read with the first column as row names.
read_counts_matrix <- function(path, delim) {
    counts <- read.table(
        path,
        sep = delim,
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )
    as.matrix(counts)
}

# Infer normal cells from the shared annotation file. Column 1 is expected to be
# the cell ID and column 2 the annotation group; groups matching ref_group_names
# are passed to CopyKAT as known normal cells.
normal_cells_from_annotations <- function(path, delim, ref_group_names) {
    if (length(ref_group_names) == 0) {
        return(character())
    }

    annotations <- read.table(
        path,
        sep = delim,
        header = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
    )

    if (ncol(annotations) < 2) {
        stop("Expected annotations_file to have at least two columns: cell_id and group.", call. = FALSE)
    }

    annotations[[1]][annotations[[2]] %in% ref_group_names]
}

# Local debug template. To run this script outside Nextflow, uncomment the block
# below, edit the paths/parameters, then run:
# Rscript modules/copycat/usr/bin/run_copycat.R
debug_args <- NULL
# debug_args <- c(
#     "--raw-counts-matrix", "/absolute/path/to/raw_counts_matrix.tsv",
#     "--annotations-file", "/absolute/path/to/cell_annotations.tsv",
#     "--out-dir", "debug_copykat",
#     "--sample-id", "copykat_debug",
#     "--ref-group-names", "normal_a,normal_b",
#     "--norm-cell-names", "",
#     "--annotations-delim", "\t",
#     "--raw-counts-delim", "\t",
#     "--id-type", "S",
#     "--cell-line", "no",
#     "--ngene-chr", "5",
#     "--low-dr", "0.05",
#     "--up-dr", "0.1",
#     "--win-size", "25",
#     "--ks-cut", "0.1",
#     "--distance", "euclidean",
#     "--output-seg", "false",
#     "--plot-genes", "true",
#     "--genome", "hg20",
#     "--num-threads", "4"
# )

# Read and normalize all process arguments before doing any heavy work.
args <- parse_args(debug_args %||% commandArgs(trailingOnly = TRUE))

# Required files and output location.
raw_counts_matrix <- required_arg(args, "raw-counts-matrix")
annotations_file <- required_arg(args, "annotations-file")
out_dir <- required_arg(args, "out-dir")
sample_id <- args[["sample-id"]] %||% tools::file_path_sans_ext(basename(raw_counts_matrix))

# Shared annotation options. Explicit normal cell names take priority over
# normals inferred from ref_group_names.
ref_group_names <- as_csv_vector(args[["ref-group-names"]])
normal_cells <- as_csv_vector(args[["norm-cell-names"]])
annotations_delim <- args[["annotations-delim"]] %||% "\t"
raw_counts_delim <- args[["raw-counts-delim"]] %||% "\t"

# CopyKAT-specific settings. Defaults mirror nextflow.config so the script can
# also be run directly for debugging.
id_type <- args[["id-type"]] %||% "S"
cell_line <- args[["cell-line"]] %||% "no"
ngene_chr <- as.integer(args[["ngene-chr"]] %||% "5")
low_dr <- as.numeric(args[["low-dr"]] %||% "0.05")
up_dr <- as.numeric(args[["up-dr"]] %||% "0.1")
win_size <- as.integer(args[["win-size"]] %||% "25")
ks_cut <- as.numeric(args[["ks-cut"]] %||% "0.1")
distance <- args[["distance"]] %||% "euclidean"
output_seg <- copykat_bool(args[["output-seg"]] %||% "false")
plot_genes <- copykat_bool(args[["plot-genes"]] %||% "true")
genome <- args[["genome"]] %||% "hg20"
num_threads <- as.integer(args[["num-threads"]] %||% "1")

load_copykat()

# If no explicit normal cells were supplied, derive them from annotation groups.
if (length(normal_cells) == 0) {
    normal_cells <- normal_cells_from_annotations(annotations_file, annotations_delim, ref_group_names)
}

# CopyKAT uses an empty string when no known normal cells are available.
norm_cell_names <- if (length(normal_cells) == 0) "" else normal_cells
rawmat <- read_counts_matrix(raw_counts_matrix, raw_counts_delim)

# CopyKAT writes several auxiliary files to the working directory, so run it
# inside the module output directory.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
old_wd <- setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

# Run CopyKAT and keep the full returned object for reproducibility.
copykat_result <- copykat(
    rawmat = rawmat,
    id.type = id_type,
    cell.line = cell_line,
    ngene.chr = ngene_chr,
    LOW.DR = low_dr,
    UP.DR = up_dr,
    win.size = win_size,
    norm.cell.names = norm_cell_names,
    KS.cut = ks_cut,
    sam.name = sample_id,
    distance = distance,
    output.seg = output_seg,
    plot.genes = plot_genes,
    genome = genome,
    n.cores = num_threads
)

# Always save the full object; write common tabular outputs when CopyKAT returns
# them for this dataset.
saveRDS(copykat_result, file = paste0(sample_id, ".copykat.rds"))

if (!is.null(copykat_result$prediction)) {
    write.table(
        copykat_result$prediction,
        file = paste0(sample_id, ".copykat.prediction.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

if (!is.null(copykat_result$CNAmat)) {
    write.table(
        copykat_result$CNAmat,
        file = paste0(sample_id, ".copykat.CNAmat.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}
