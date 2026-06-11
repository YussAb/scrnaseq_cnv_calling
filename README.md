# Minimal inferCNV Nextflow Pipeline

This is a compact DSL2 pipeline with one local module:

```text
.
├── main.nf
├── nextflow.config
└── modules
    └── infercnv
        ├── main.nf
        └── usr
            └── bin
                └── run_infercnv.R
```

The structure borrows the useful, lightweight parts of nf-core style: parameters are centralised in `nextflow.config`, the workflow includes a local module, inputs are passed with a `meta` map, the process emits `versions.yml`, and runtime reporting is enabled.

## Run

```bash
nextflow run . \
    --raw_counts_matrix path/to/raw_counts_matrix.tsv \
    --annotations_file path/to/cell_annotations.tsv \
    --gene_order_file path/to/gene_order.tsv \
    --ref_group_names normal_cells \
    --outdir results
```

Use a comma-separated string or a list in config for multiple reference groups:

```bash
nextflow run . \
    --raw_counts_matrix path/to/raw_counts_matrix.tsv \
    --annotations_file path/to/cell_annotations.tsv \
    --gene_order_file path/to/gene_order.tsv \
    --ref_group_names normal_a,normal_b
```

The execution environment must provide R and the `infercnv` R package. You can set `params.infercnv_container` in `nextflow.config` or at runtime and enable a container profile:

```bash
nextflow run . -profile docker \
    --infercnv_container '<image-with-r-and-infercnv>' \
    --raw_counts_matrix path/to/raw_counts_matrix.tsv \
    --annotations_file path/to/cell_annotations.tsv \
    --gene_order_file path/to/gene_order.tsv
```
