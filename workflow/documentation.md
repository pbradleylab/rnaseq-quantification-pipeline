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
- `rules/metrics.smk`: runs FastQC, FastQ Screen, GUNC, strandedness checks,
  MultiQC, the final per-sample QC summary, and the sample identity report.
- `rules/quantification.smk`: runs Kallisto, STAR, or FeatureCounts
  quantification and merges per-sample outputs.
- `rules/deseq.smk`: runs DESeq2 and writes differential-expression tables plus
  independent plot outputs under plot-specific results directories.

## Script Helpers

- `scripts/utils.py`: shared helper functions for PEP metadata, read routing,
  and project-level rule inputs.
- `scripts/combine_star_counts.py`: merges STAR `ReadsPerGene.out.tab` files.
  The `--strandedness` option selects the STAR count column.
- `scripts/check_strandedness.py`: infers library strandedness from STAR
  `ReadsPerGene.out.tab` columns and compares it to the configured STAR or
  FeatureCounts strandedness setting.
- `scripts/merge_strandedness_reports.py`: merges per-sample strandedness
  reports into the final project-level QC table.
- `scripts/combine_featurecounts.py`: merges per-sample FeatureCounts outputs.
- `scripts/combine_kallisto.py`: merges Kallisto `abundance.tsv` files from
  paired-end and single-end output directories.
- `scripts/summarize_qc.py`: builds the final per-sample QC table.
- `scripts/sample_identity_report.py`: fingerprints expression profiles from
  the merged count matrix, reports nearest expression neighbors, and flags
  low-similarity or metadata-discordant samples.
- `scripts/check_count_annotation_overlap.py`: compares count-matrix
  `target_id` values with configured annotation feature IDs before DESeq2 runs.
- `scripts/summarize_gene_biotype_counts.py`: summarizes count support by
  `gene_biotype` and count class when the annotation contains `gene_biotype`.
- `scripts/DESeq.R`: runs DESeq2 in a selected mode to generate the result
  table, significant gene tables, normalized count matrix, transformed count
  matrix, Cook's distance reports, volcano plot, MA plot, normalized expression
  plots, sample distance heatmap, PCA plot, or library size and size-factor
  plot.

## Important Outputs

- `results/{project}/final/multiqc/multiqc_report.html`
- `results/{project}/final/qc/{project}_sample_qc_summary.tsv`
- `results/{project}/final/qc/{project}_sample_identity_report.tsv`
- `results/{project}/final/qc/{project}_strandedness_summary.tsv`
- `results/{project}/count_annotation_overlap/{project}_count_annotation_overlap.tsv`
- `results/{project}/gene_biotype_count_summary/{project}_gene_biotype_count_summary.tsv`
- `results/{project}/final/quantification/star/gene_counts_all_samples.tsv`
- `results/{project}/final/quantification/featurecounts/gene_counts_all_samples.tsv`
- `results/{project}/final/quantification/kallisto/transcript_counts_all_samples.tsv`
- `results/{project}/deseq2/{project}.tsv`
- `results/{project}/significant_gene_tables/{project}_significant_genes.tsv`
- `results/{project}/significant_gene_tables/{project}_significant_upregulated_genes.tsv`
- `results/{project}/significant_gene_tables/{project}_significant_downregulated_genes.tsv`
- `results/{project}/normalized_counts/{project}_normalized_counts.tsv`
- `results/{project}/transformed_counts/{project}_transformed_counts.tsv`
- `results/{project}/cooks_reports/{project}_cooks_gene_report.tsv`
- `results/{project}/cooks_reports/{project}_cooks_sample_report.tsv`
- `results/{project}/volcano_plot/{project}_volcano_plot.png`
- `results/{project}/volcano_plot/{project}_volcano_plot.svg`
- `results/{project}/ma_plot/{project}_ma_plot.png`
- `results/{project}/ma_plot/{project}_ma_plot.svg`
- `results/{project}/sample_distance_heatmap/{project}_sample_distance_heatmap.png`
- `results/{project}/sample_distance_heatmap/{project}_sample_distance_heatmap.svg`
- `results/{project}/pca/{project}_pca.png`
- `results/{project}/pca/{project}_pca.svg`
- `results/{project}/library_sizes_size_factors/{project}_library_sizes_size_factors.png`
- `results/{project}/library_sizes_size_factors/{project}_library_sizes_size_factors.svg`

## Strandedness Checks

`strandedness_check.enabled` controls whether STAR/FeatureCounts workflows infer
library strandedness from STAR `ReadsPerGene.out.tab` outputs. The checker
compares the inferred value to `star.gene_counts_strandedness` for STAR
quantification or `featurecounts.strandedness` for FeatureCounts quantification.
Set `strandedness_check.mode` to `warn` to write mismatch reports without
failing the workflow, or `fail` to stop when inferred strandedness conflicts
with the configured count mode.

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
