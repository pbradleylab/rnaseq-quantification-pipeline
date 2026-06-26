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
- `rules/deseq.smk`: runs DESeq2 and writes differential-expression tables plus
  independent plot outputs under the differential-expression results directory.

## Script Helpers

- `scripts/utils.py`: shared helper functions for PEP metadata, read routing,
  and project-level rule inputs.
- `scripts/combine_star_counts.py`: merges STAR `ReadsPerGene.out.tab` files.
  The `--strandedness` option selects the STAR count column.
- `scripts/combine_featurecounts.py`: merges per-sample FeatureCounts outputs.
- `scripts/combine_kallisto.py`: merges Kallisto `abundance.tsv` files from
  paired-end and single-end output directories.
- `scripts/summarize_qc.py`: builds the final per-sample QC table.
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
- `results/{project}/final/qc/{project}_count_annotation_overlap.tsv`
- `results/{project}/final/qc/{project}_gene_biotype_count_summary.tsv`
- `results/{project}/final/quantification/star/gene_counts_all_samples.tsv`
- `results/{project}/final/quantification/featurecounts/gene_counts_all_samples.tsv`
- `results/{project}/final/quantification/kallisto/transcript_counts_all_samples.tsv`
- `results/{project}/differential_expression/{project}.tsv`
- `results/{project}/differential_expression/{project}_significant_genes.tsv`
- `results/{project}/differential_expression/{project}_significant_upregulated_genes.tsv`
- `results/{project}/differential_expression/{project}_significant_downregulated_genes.tsv`
- `results/{project}/differential_expression/{project}_normalized_counts.tsv`
- `results/{project}/differential_expression/{project}_transformed_counts.tsv`
- `results/{project}/differential_expression/{project}_cooks_gene_report.tsv`
- `results/{project}/differential_expression/{project}_cooks_sample_report.tsv`
- `results/{project}/differential_expression/{project}_volcano_plot.png`
- `results/{project}/differential_expression/{project}_volcano_plot.svg`
- `results/{project}/differential_expression/{project}_ma_plot.png`
- `results/{project}/differential_expression/{project}_ma_plot.svg`
- `results/{project}/differential_expression/{project}_sample_distance_heatmap.png`
- `results/{project}/differential_expression/{project}_sample_distance_heatmap.svg`
- `results/{project}/differential_expression/{project}_pca.png`
- `results/{project}/differential_expression/{project}_pca.svg`
- `results/{project}/differential_expression/{project}_library_sizes_size_factors.png`
- `results/{project}/differential_expression/{project}_library_sizes_size_factors.svg`

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
