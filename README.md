# Minimal inferCNV Nextflow Pipeline

This is a compact DSL2 pipeline with two local CNV modules:

```text
.
├── main.nf
├── nextflow.config
└── modules
    ├── copycat
    │   ├── main.nf
    │   └── usr
    │       └── bin
    │           └── run_copycat.R
    └── infercnv
        ├── main.nf
        └── usr
            └── bin
                └── run_infercnv.R
```

The structure borrows the useful, lightweight parts of nf-core style: parameters are centralised in `nextflow.config`, the workflow includes local modules, inputs are passed with a `meta` map, each process emits `versions.yml`, and runtime reporting is enabled.

The CopyCAT module wraps the R package `copykat`. It reuses the same raw counts matrix and annotations file used by inferCNV. When `ref_group_names` is set, matching cells from `annotations_file` are passed to CopyKAT as known normal cells through `norm.cell.names`.

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

Run only one module when needed:

```bash
nextflow run . -params-file params.yaml --run_copycat false
nextflow run . -params-file params.yaml --run_infercnv false
```
