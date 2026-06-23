# RNA-seq Quantification Workflow

This Snakemake workflow processes bulk RNA-seq data from FASTQ files through
trimming, QC, quantification, MultiQC reporting, and basic DESeq2 differential
expression outputs. It supports mixed paired-end and single-end projects.

## Main Outputs

For each project, the workflow writes final outputs under `results/{project}/`.

- `final/multiqc/multiqc_report.html`: combined QC report.
- `final/qc/{project}_sample_qc_summary.tsv`: per-sample QC summary table.
- `final/quantification/star/gene_counts_all_samples.tsv`: STAR gene-count matrix when `quantification_tool` is `star`.
- `final/quantification/kallisto/transcript_counts_all_samples.tsv`: Kallisto count matrix when `quantification_tool` is `kallisto`.
- `final/quantification/kallisto/transcript_tpms_all_samples.tsv`: Kallisto TPM matrix.
- `differential_expression/{project}.tsv`: DESeq2 results table.
- `differential_expression/{project}.png`, `.svg`, `.pdf`: volcano plot.
- `differential_expression/{project}_normalized_expression_boxplot.*`: normalized expression boxplot in PNG, SVG, and PDF.
- `differential_expression/{project}_normalized_expression_density.*`: normalized expression density plot in PNG, SVG, and PDF.

## Installation

Install conda or mamba first. Then create an environment with Snakemake:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

The workflow-specific tool environments are defined in `workflow/envs/` and are
created automatically when running Snakemake with `--use-conda`.

## Test Data

The repository includes a small test reference and test sample sheets. To check
that the workflow DAG builds:

```bash
conda activate snakemake
snakemake --use-conda --cores 1 -n
```

The test FASTQs can be generated from the `test/` directory:

```bash
cd test
snakemake --use-conda --cores 4
cd ..
```

## Input Configuration

The main configuration files are:

- `config/config.json`: project-level settings, reference paths, metadata, and quantification choice.
- `config/tools.json`: tool resource settings and analysis parameters.
- `config/pepconfig.yml`: PEP sample table configuration.
- `test/test_samples.csv` and `test/test_subsamples.csv`: example sample sheets.

### Sample Table

The sample table must include:

- `sample_name`: sample identifier. This must match the sample IDs in metadata and the subsample table.
- `alternate_id`: optional alternate identifier.
- `project`: project name used in output paths.
- `organism`: organism or strain.
- `sample_type`: use `rna`.
- `seq_method`: `paired_end` or `single_end`.

### Subsample Table

The subsample table must include:

- `sample_name`: sample identifier from the sample table.
- `subsample`: subsample identifier used to discover FASTQ files.
- `protocol`: use `rna`.
- `seq_method`: `paired_end` or `single_end`.

FASTQ discovery is configured in `config/pepconfig.yml`. By default, reads are
derived from:

```yaml
sample_modifiers:
  append:
    raw_data: "test/Mycolicibacterium_smegmatis/"
  derive:
    attributes: [reads]
    sources:
      rna: "{raw_data}/{subsample}*.fastq.gz"
```

For paired-end data, read names should contain a recognizable R1/R2 marker such
as `_R1`, `_R2`, `.R1`, `.R2`, `_1.fastq`, or `_2.fastq`.

### Reference Configuration

In `config/config.json`, set:

- `genome.is_local`: `true` when using a local FASTA.
- `genome.name`: stable reference name used in output paths.
- `genome.path`: local FASTA path when `is_local` is `true`.
- `genome.url`: remote FASTA URL when `is_local` is `false`.
- `gff3.is_local`: `true` when using a local annotation file.
- `gff3.path`: local GFF/GTF path when `is_local` is `true`.
- `gff3.url`: remote annotation URL when `is_local` is `false`.
- `metadata`: DESeq2 metadata CSV.
- `run_gunc`: `"true"` or `"false"`.
- `quantification_tool`: `"star"` or `"kallisto"`.

## Tool Configuration

Thread and memory settings are in `config/tools.json`.

Important STAR options:

- `star.read_files`: command used by STAR to read compressed FASTQ files, for example `gunzip -c`.
- `star.type`: STAR BAM output type, for example `BAM Unsorted`.
- `star.gene_counts_strandedness`: count column used from `ReadsPerGene.out.tab`.
  - `unstranded`: STAR column 2.
  - `forward`: STAR column 3.
  - `reverse`: STAR column 4.
- `star.other`: extra STAR arguments.

DESeq2 options:

- `deseq2.variable_to_analyze`: metadata column used in the design formula.
- `deseq2.reference_in_variable`: reference level for that metadata column.
- `deseq2.log2fc_threshold`: volcano plot fold-change threshold.
- `deseq2.padj_threshold`: adjusted p-value threshold.
- `deseq2.label_top_n`: number of top genes to label on the volcano plot.

## Running The Workflow

From the repository root:

```bash
conda activate snakemake
snakemake --use-conda --cores 4
```

For a dry run:

```bash
snakemake --use-conda --cores 1 -n
```

On SLURM with Snakemake's executor support, adapt this pattern to your cluster:

```bash
snakemake --use-conda --slurm --default-resources slurm_account=YOUR_ACCOUNT --jobs 50
```

## QC And Reports

The workflow runs:

- Trim Galore for read trimming.
- FastQC on raw and trimmed reads.
- FastQ Screen for contamination screening.
- STAR or Kallisto for quantification.
- GUNC on the reference genome when enabled.
- MultiQC for combined reporting.
- DESeq2 for differential expression and plots.

The per-sample QC summary combines key metadata, FastQC status, FastQ Screen
report presence, assigned counts, and STAR mapping metrics when STAR is used.

## Common Problems

- **Metadata/count mismatch in DESeq2**: sample IDs in the metadata must match
  count-matrix columns exactly. The workflow now stops with a clear error when
  they do not match.
- **Paired-end FASTQs not found**: check R1/R2 naming. The workflow recognizes
  common endings such as `_R1`, `_R2`, `_1.fastq`, and `_2.fastq`.
- **Wrong STAR count column**: set `star.gene_counts_strandedness` to
  `unstranded`, `forward`, or `reverse` to match the library strandedness.
- **Run location**: run the main workflow from the repository root, not from
  `test/` or `workflow/`.

## DAG

An example workflow DAG is available at:

![Workflow DAG](images/dag.svg)
