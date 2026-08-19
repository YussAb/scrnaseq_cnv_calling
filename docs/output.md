# Output

Runtime reports are written under:

```text
<outdir>/pipeline_info/
```

inferCNV results are published under:

```text
<outdir>/infercnv/infercnv/
```

CopyKAT results are published under:

```text
<outdir>/copykat/copykat/
```

Avoid setting `outdir` to a path that already ends with the module name unless
that extra nesting is intentional.
