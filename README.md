# Minimal inferCNV / CopyKAT Nextflow Pipeline

This is a compact DSL2 pipeline for running CNV inference with two local modules:

- `inferCNV`, through the R package `infercnv`
- `CopyKAT`, through the R package `copykat`

The workflow reuses the same raw counts matrix and cell annotation file for both
methods. When `ref_group_names` is provided, the same reference groups are used by
inferCNV and are also converted into known normal cell names for CopyKAT unless
`copykat_norm_cell_names` is set explicitly.

```text
.
├── main.nf
├── nextflow.config
├── infercnv.params.yaml
├── copycat.params.yaml
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

## Run

Run both modules:

```bash
nextflow run . \
    --raw_counts_matrix path/to/raw_counts_matrix.tsv \
    --annotations_file path/to/cell_annotations.tsv \
    --gene_order_file path/to/gene_order.tsv \
    --ref_group_names normal_cells \
    --outdir results
```

Use comma-separated reference groups when more than one annotation group should be
treated as normal/reference cells:

```bash
nextflow run . \
    --raw_counts_matrix path/to/raw_counts_matrix.tsv \
    --annotations_file path/to/cell_annotations.tsv \
    --gene_order_file path/to/gene_order.tsv \
    --ref_group_names normal_a,normal_b
```

Run only inferCNV:

```bash
nextflow run . -params-file infercnv.params.yaml
```

Run only CopyKAT:

```bash
nextflow run . -params-file copycat.params.yaml
```

You can also disable either module from the command line:

```bash
nextflow run . -params-file params.yaml --run_copykat false
nextflow run . -params-file params.yaml --run_infercnv false
```

The pipeline defines Docker and Singularity profiles. The default containers are
configured in `nextflow.config`.

```bash
nextflow run . -profile docker -params-file infercnv.params.yaml
nextflow run . -profile docker -params-file copycat.params.yaml
```

## Inputs

| Input | Meaning |
| --- | --- |
| `raw_counts_matrix` | Raw gene-by-cell count matrix. The first column is used as gene names and the remaining columns are cell counts. |
| `annotations_file` | Cell annotation file used by both tools. It should contain at least two columns: cell ID and group name. |
| `gene_order_file` | Gene genomic position file required by inferCNV. This is only required when `run_infercnv` is true. |

## Shared Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `run_infercnv` | `true` | Enables the inferCNV module. At least one of `run_infercnv` or `run_copykat` must be true. |
| `run_copykat` | `true` | Enables the CopyKAT module. |
| `raw_counts_matrix` | `null` | Path to the raw counts matrix. Required for both modules. |
| `annotations_file` | `null` | Path to the cell annotation file. Required for both modules. |
| `ref_group_names` | `[]` | Reference or normal annotation group name(s). Use a comma-separated string such as `normal_a,normal_b` or a list in a params file. |
| `sample_id` | `null` | Sample name used by CopyKAT output files. If unset, the raw counts matrix basename is used. |
| `annotations_delim` | `'\t'` | Delimiter used to read `annotations_file`. |
| `outdir` | `results` | Output directory. |
| `publish_dir_mode` | `copy` | Nextflow publish mode for result files. |

## inferCNV Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `gene_order_file` | `null` | Gene order file passed to `CreateInfercnvObject`. Required when inferCNV is enabled. |
| `cutoff` | `0.1` | Minimum average read count used by `infercnv::run`. A common value is `0.1` for droplet-based single-cell data. |
| `cluster_by_groups` | `true` | Tells inferCNV to cluster cells by the annotation groups. |
| `denoise` | `true` | Enables inferCNV denoising. |
| `hmm` | `false` | Enables the inferCNV HMM step when true. |

inferCNV outputs are published under:

```text
results/infercnv
```

## CopyKAT Parameters

The CopyKAT module uses `raw_counts_matrix` and `annotations_file`. If
`copykat_norm_cell_names` is empty, the wrapper reads `annotations_file` and uses
cells whose annotation group is listed in `ref_group_names` as known normal cells.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `copykat_raw_counts_delim` | `'\t'` | Delimiter used to read `raw_counts_matrix` for CopyKAT. |
| `copykat_norm_cell_names` | `[]` | Explicit normal cell IDs passed to CopyKAT as `norm.cell.names`. When set, this overrides normals inferred from `ref_group_names`. |
| `copykat_id_type` | `S` | Gene ID type passed as `id.type`. `S` is commonly used for gene symbols. |
| `copykat_cell_line` | `no` | Passed as `cell.line`; set according to whether the input is from a cell line experiment. |
| `copykat_ngene_chr` | `5` | Minimum number of genes required per chromosome, passed as `ngene.chr`. |
| `copykat_low_dr` | `0.05` | Lower dropout-rate filter passed as `LOW.DR`. |
| `copykat_up_dr` | `0.1` | Upper dropout-rate filter passed as `UP.DR`. |
| `copykat_win_size` | `25` | Genomic window size passed as `win.size`. |
| `copykat_ks_cut` | `0.1` | Segmentation threshold passed as `KS.cut`. |
| `copykat_distance` | `euclidean` | Distance metric passed to CopyKAT clustering. |
| `copykat_output_seg` | `false` | Writes CopyKAT segmentation output when true. |
| `copykat_plot_genes` | `true` | Enables gene-level plotting in CopyKAT. |
| `copykat_genome` | `hg20` | Genome build passed to CopyKAT. |

CopyKAT outputs are published under:

```text
results/copykat
```

The wrapper always saves the full CopyKAT result object as:

```text
<sample_id>.copykat.rds
```

When available, it also writes:

```text
<sample_id>.copykat.prediction.tsv
<sample_id>.copykat.CNAmat.tsv
```

## Runtime Reports

Nextflow runtime reports are written to:

```text
results/pipeline_info/
```

This includes `timeline.html`, `report.html`, `trace.txt`, and
`pipeline_dag.html`.
