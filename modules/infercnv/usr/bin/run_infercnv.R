#!/usr/bin/env Rscript

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

required_arg <- function(args, name) {
    value <- args[[name]]
    if (is.null(value) || identical(value, "")) {
        stop("Missing required argument --", name, call. = FALSE)
    }
    value
}

as_bool <- function(value) {
    tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

`%||%` <- function(lhs, rhs) {
    if (is.null(lhs)) rhs else lhs
}

as_ref_groups <- function(value) {
    if (is.null(value) || isTRUE(value) || identical(value, "")) {
        return(NULL)
    }
    trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

raw_counts_matrix <- required_arg(args, "raw-counts-matrix")
annotations_file <- required_arg(args, "annotations-file")
gene_order_file <- required_arg(args, "gene-order-file")
out_dir <- required_arg(args, "out-dir")

annotations_delim <- args[["annotations-delim"]] %||% "\t"
ref_group_names <- as_ref_groups(args[["ref-group-names"]])
cutoff <- as.numeric(args[["cutoff"]] %||% "0.1")
cluster_by_groups <- as_bool(args[["cluster-by-groups"]] %||% "true")
denoise <- as_bool(args[["denoise"]] %||% "true")
hmm <- as_bool(args[["hmm"]] %||% "false")
num_threads <- as.integer(args[["num-threads"]] %||% "1")

if (!requireNamespace("infercnv", quietly = TRUE)) {
    stop(
        "The R package 'infercnv' is not installed. ",
        "Run with a container/conda environment that provides infercnv.",
        call. = FALSE
    )
}

infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = raw_counts_matrix,
    annotations_file = annotations_file,
    delim = annotations_delim,
    gene_order_file = gene_order_file,
    ref_group_names = ref_group_names
)

infercnv::run(
    infercnv_obj = infercnv_obj,
    cutoff = cutoff,
    out_dir = out_dir,
    cluster_by_groups = cluster_by_groups,
    denoise = denoise,
    HMM = hmm,
    num_threads = num_threads
)
