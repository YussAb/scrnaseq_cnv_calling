#!/usr/bin/env Rscript

# Lightweight inferCNV runner used by the Nextflow module. It keeps the process
# script small and maps Nextflow parameters to the infercnv R API.

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

# Convert Nextflow-style boolean values to R logical values.
as_bool <- function(value) {
    tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

# Small default-value helper: lhs %||% rhs returns rhs only when lhs is NULL.
`%||%` <- function(lhs, rhs) {
    if (is.null(lhs)) rhs else lhs
}

# inferCNV accepts NULL when no reference groups are supplied; otherwise it
# expects a character vector of annotation group names.
as_ref_groups <- function(value) {
    if (is.null(value) || isTRUE(value) || identical(value, "")) {
        return(NULL)
    }
    trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

# Local debug template. To run this script outside Nextflow, uncomment the block
# below, edit the paths/parameters, then run:
# Rscript modules/infercnv/usr/bin/run_infercnv.R
debug_args <- NULL
# debug_args <- c(
#     "--raw-counts-matrix", "/absolute/path/to/raw_counts_matrix.tsv",
#     "--annotations-file", "/absolute/path/to/cell_annotations.tsv",
#     "--gene-order-file", "/absolute/path/to/gene_order.tsv",
#     "--out-dir", "debug_infercnv",
#     "--annotations-delim", "\t",
#     "--ref-group-names", "normal_a,normal_b",
#     "--cutoff", "0.1",
#     "--cluster-by-groups", "true",
#     "--denoise", "true",
#     "--hmm", "false",
#     "--num-threads", "4"
# )

# Read and normalize all process arguments before creating the inferCNV object.
args <- parse_args(debug_args %||% commandArgs(trailingOnly = TRUE))

# Required input files and output directory.
raw_counts_matrix <- required_arg(args, "raw-counts-matrix")
annotations_file <- required_arg(args, "annotations-file")
gene_order_file <- required_arg(args, "gene-order-file")
out_dir <- required_arg(args, "out-dir")

# inferCNV options. Defaults mirror nextflow.config so the script can also be
# run directly for debugging.
annotations_delim <- args[["annotations-delim"]] %||% "\t"
ref_group_names <- as_ref_groups(args[["ref-group-names"]])
cutoff <- as.numeric(args[["cutoff"]] %||% "0.1")
cluster_by_groups <- as_bool(args[["cluster-by-groups"]] %||% "true")
denoise <- as_bool(args[["denoise"]] %||% "true")
hmm <- as_bool(args[["hmm"]] %||% "false")
num_threads <- as.integer(args[["num-threads"]] %||% "1")

# Fail early with a clear message if the selected container/environment is wrong.
if (!requireNamespace("infercnv", quietly = TRUE)) {
    stop(
        "The R package 'infercnv' is not installed. ",
        "Run with a container/conda environment that provides infercnv.",
        call. = FALSE
    )
}

# Create the inferCNV object from the shared pipeline inputs. The annotation file
# maps cells to groups, and the gene order file maps genes to genomic positions.
infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = raw_counts_matrix,
    annotations_file = annotations_file,
    delim = annotations_delim,
    gene_order_file = gene_order_file,
    ref_group_names = ref_group_names
)

# Run the main inferCNV workflow with the requested denoising, clustering, and
# optional HMM settings.
infercnv::run(
    infercnv_obj = infercnv_obj,
    cutoff = cutoff,
    out_dir = out_dir,
    cluster_by_groups = cluster_by_groups,
    denoise = denoise,
    HMM = hmm,
    num_threads = num_threads
)
