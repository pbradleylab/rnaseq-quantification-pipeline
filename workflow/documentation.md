# Workflow Documentation

This directory contains the Snakemake implementation for the RNA-seq
quantification workflow. Run the workflow from the repository root, not from
inside `workflow/`.

## Rule Groups

- `rules/utilities.smk`: symlinks local references and input FASTQs into stable
  workflow paths.
- `rules/resources.smk`: prepares reference resources, Kallisto indices, STAR
  indices, and optional GUNC databases.
- `rules/trimming.smk`: trims paired-end and single-end reads with Trim Galore.
- `rules/metrics.smk`: runs FastQC, FastQ Screen, GUNC, MultiQC, and the final
  per-sample QC summary.
- `rules/quantification.smk`: runs Kallisto, STAR, or FeatureCounts
  quantification and merges per-sample outputs.
- `rules/deseq.smk`: runs DESeq2 and writes differential-expression tables and
  plots.

## Script Helpers

- `scripts/utils.py`: shared helper functions for PEP metadata, read routing,
  and project-level rule inputs.
- `scripts/combine_star_counts.py`: merges STAR `ReadsPerGene.out.tab` files.
  The `--strandedness` option selects the STAR count column.
- `scripts/combine_featurecounts.py`: merges per-sample FeatureCounts outputs.
- `scripts/combine_kallisto.py`: merges Kallisto `abundance.tsv` files from
  paired-end and single-end output directories.
- `scripts/summarize_qc.py`: builds the final per-sample QC table.
- `scripts/DESeq.R`: runs DESeq2 and generates result tables plus PNG, SVG, and
  PDF plots.

## Important Outputs

- `results/{project}/final/multiqc/multiqc_report.html`
- `results/{project}/final/qc/{project}_sample_qc_summary.tsv`
- `results/{project}/final/quantification/star/gene_counts_all_samples.tsv`
- `results/{project}/final/quantification/featurecounts/gene_counts_all_samples.tsv`
- `results/{project}/final/quantification/kallisto/transcript_counts_all_samples.tsv`
- `results/{project}/differential_expression/{project}.tsv`

## Development Checks

Use the Snakemake conda environment from the repository root:

```bash
conda activate snakemake
snakemake --use-conda --cores 1 -n
snakemake --lint
```

The lint command can report style warnings for hardcoded output directories and
helper functions inside rule files. Treat DAG failures and script errors as
higher priority than those style warnings.
