# Usage

Run the pipeline with either a parameter file or explicit command-line options.

```bash
nextflow run . \
    --raw_counts_matrix path/to/raw_counts_matrix.tsv \
    --annotations_file path/to/cell_annotations.tsv \
    --gene_order_file path/to/gene_order.tsv \
    --ref_group_names normal_cells \
    --outdir results
```

Run only inferCNV:

```bash
nextflow run . -params-file infercnv.params.yaml
```

Run only CopyKAT:

```bash
nextflow run . -params-file copykat.params.yaml
```

The included parameter files are local examples and may contain absolute paths
for this environment. Update input paths and `outdir` before running elsewhere.
