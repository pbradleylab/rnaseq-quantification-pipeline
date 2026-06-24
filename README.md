# RNA-seq Quantification Workflow

This repository contains a Snakemake workflow for bulk RNA-seq processing,
quality control, quantification, reporting, and basic DESeq2 differential
expression analysis. It was originally built around isolate RNA-seq workflows,
but now supports STAR, FeatureCounts, and Kallisto quantification modes.

For MAG-derived samples or projects with limited annotations, an alignment-based
approach is usually a better first choice than transcriptome-only
quantification. If the annotation set is sparse, consider improving annotation
before using this workflow.

The workflow supports projects that contain paired-end reads, single-end reads,
or a mixture of both.

## What The Workflow Can Do

### Input Handling

- Read sample and subsample metadata through a PEP configuration.
- Track sample-level metadata, project names, organisms, sequencing method, and
  subsample identifiers.
- Handle paired-end and single-end RNA-seq libraries in the same project.
- Discover FASTQ files from configurable path patterns in `config/pepconfig.yml`.
- Recognize common paired-end read name patterns such as `_R1`, `_R2`,
  `.R1`, `.R2`, `_1.fastq`, and `_2.fastq`.
- Symlink raw FASTQs into stable workflow-managed paths under `resources/`.

### Reference Preparation

- Use local genome FASTA and annotation files.
- Download remote genome FASTA and annotation files when configured.
- Symlink local references into reproducible workflow paths.
- Convert GFF/GFF3 annotations to GTF for STAR and FeatureCounts-related steps.
- Generate a transcript FASTA from the genome and annotation with `gffread`.
- Build a STAR genome index for alignment-based quantification.
- Build a Kallisto transcriptome index for transcript quantification.
- Optionally download the GUNC database for reference quality checks.

### Read Processing

- Trim paired-end reads with Trim Galore.
- Trim single-end reads with Trim Galore.
- Run FastQC on raw reads.
- Run FastQC on trimmed reads.
- Run FastQ Screen to screen reads against contaminant/reference databases.
- Include raw QC, trimmed QC, contamination screening, alignment outputs, and
  quantification outputs in MultiQC.

### Quantification

The quantification mode is controlled by `quantification_tool` in
`config/config.json`.

Supported values are:

- `star`
- `featurecounts`
- `kallisto`

#### STAR Gene Counts

When `quantification_tool` is set to `star`, the workflow can:

- Align paired-end and single-end reads with STAR.
- Produce STAR gene count files with `--quantMode GeneCounts`.
- Produce STAR transcriptome BAMs with `--quantMode TranscriptomeSAM`.
- Merge per-sample STAR gene counts into
  `results/{project}/final/quantification/star/gene_counts_all_samples.tsv`.
- Select the STAR count column by library strandedness using
  `star.gene_counts_strandedness` in `config/tools.json`.

Supported STAR count modes are:

- `unstranded`: use STAR column 2.
- `forward`: use STAR column 3.
- `reverse`: use STAR column 4.

#### FeatureCounts Gene Counts

When `quantification_tool` is set to `featurecounts`, the workflow can:

- Align paired-end and single-end reads with STAR.
- Run FeatureCounts on the resulting BAM files.
- Count configurable annotation features such as `exon` or `CDS`.
- Use a configurable annotation attribute such as `gene_id`.
- Support unstranded, stranded, and reversely stranded counting through the
  FeatureCounts `-s` setting.
- Merge per-sample FeatureCounts outputs into
  `results/{project}/final/quantification/featurecounts/gene_counts_all_samples.tsv`.
- Include FeatureCounts summary files in MultiQC.

#### Kallisto Transcript Quantification

When `quantification_tool` is set to `kallisto`, the workflow can:

- Build a Kallisto transcriptome index from the generated transcript FASTA.
- Quantify paired-end reads with Kallisto.
- Quantify single-end reads with Kallisto using configured single-end estimates.
- Merge transcript estimated counts into
  `results/{project}/final/quantification/kallisto/transcript_counts_all_samples.tsv`.
- Merge transcript TPMs into
  `results/{project}/final/quantification/kallisto/transcript_tpms_all_samples.tsv`.

### Quality Reports

For each project, the workflow produces:

- A MultiQC HTML report.
- A per-sample QC summary table.
- FastQC reports for raw reads.
- FastQC reports for trimmed reads.
- FastQ Screen reports.
- STAR mapping logs when STAR is used directly or as the FeatureCounts aligner.
- FeatureCounts assignment summaries when FeatureCounts is used.
- Optional GUNC reports for reference genome quality checks.

The final QC summary table includes metadata columns, DESeq2 inclusion status,
assigned count totals, FastQC/FastQ Screen report presence, and STAR mapping
metrics when available.

### Differential Expression

The workflow runs a basic DESeq2 analysis after quantification. It can:

- Read the merged count matrix from STAR, FeatureCounts, or Kallisto.
- Read sample metadata from the CSV configured in `config/config.json`.
- Match count matrix columns to metadata sample IDs.
- Stop with a clear error when count and metadata sample names do not match.
- Use a configurable metadata variable as the DESeq2 design factor.
- Use a configurable reference level for that variable.
- Export a DESeq2 result table.
- Generate a volcano plot.
- Generate normalized expression boxplots.
- Generate normalized expression density plots.
- Save generated plots as PNG, SVG, and PDF.

Current DESeq2 configuration is controlled by `config/tools.json`:

- `deseq2.variable_to_analyze`
- `deseq2.reference_in_variable`
- `deseq2.log2fc_threshold`
- `deseq2.padj_threshold`
- `deseq2.label_top_n`

## Main Outputs

For each project, final outputs are written under `results/{project}/`.

- `final/multiqc/multiqc_report.html`: combined QC report.
- `final/qc/{project}_sample_qc_summary.tsv`: per-sample QC summary table.
- `final/quantification/star/gene_counts_all_samples.tsv`: STAR gene-count
  matrix when `quantification_tool` is `star`.
- `final/quantification/featurecounts/gene_counts_all_samples.tsv`:
  FeatureCounts gene-count matrix when `quantification_tool` is
  `featurecounts`.
- `final/quantification/kallisto/transcript_counts_all_samples.tsv`: Kallisto
  transcript estimated-count matrix when `quantification_tool` is `kallisto`.
- `final/quantification/kallisto/transcript_tpms_all_samples.tsv`: Kallisto
  transcript TPM matrix when `quantification_tool` is `kallisto`.
- `differential_expression/{project}.tsv`: DESeq2 results table.
- `differential_expression/{project}.png`: volcano plot.
- `differential_expression/{project}.svg`: volcano plot.
- `differential_expression/{project}.pdf`: volcano plot.
- `differential_expression/{project}_normalized_expression_boxplot.png`
- `differential_expression/{project}_normalized_expression_boxplot.svg`
- `differential_expression/{project}_normalized_expression_boxplot.pdf`
- `differential_expression/{project}_normalized_expression_density.png`
- `differential_expression/{project}_normalized_expression_density.svg`
- `differential_expression/{project}_normalized_expression_density.pdf`

Intermediate outputs are written under `resources/`, `results/{project}/`, and
`logs/`. Snakemake-managed Conda environments are written under `.snakemake/`.

## Repository Layout

- `workflow/Snakefile`: top-level Snakemake entry point.
- `workflow/rules/utilities.smk`: symlinks local references and FASTQs.
- `workflow/rules/resources.smk`: prepares reference resources and indices.
- `workflow/rules/trimming.smk`: trims paired-end and single-end reads.
- `workflow/rules/metrics.smk`: runs FastQC, FastQ Screen, GUNC, MultiQC, and
  QC summary generation.
- `workflow/rules/quantification.smk`: runs STAR, FeatureCounts, or Kallisto.
- `workflow/rules/deseq.smk`: runs DESeq2 and differential-expression plots.
- `workflow/scripts/`: helper scripts for count merging, QC summaries, and
  DESeq2 plotting.
- `workflow/envs/`: Conda environment YAML files used by Snakemake.
- `config/config.json`: project-level workflow settings.
- `config/tools.json`: tool settings, resources, and analysis parameters.
- `config/pepconfig.yml`: PEP configuration for sample discovery.
- `test/`: small example data and metadata.

## Installation

Install Conda or Mamba first. Then create and activate a Snakemake environment:

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

Then run the main workflow from the repository root:

```bash
snakemake --use-conda --cores 4
```

## Input Configuration

The main configuration files are:

- `config/config.json`: project-level settings, reference paths, metadata, and
  quantification choice.
- `config/tools.json`: tool resource settings and analysis parameters.
- `config/pepconfig.yml`: PEP sample table configuration.
- `test/test_samples.csv` and `test/test_subsamples.csv`: example sample
  sheets.

### Sample Table

The sample table must include:

- `sample_name`: sample identifier. This must match the sample IDs in metadata
  and the subsample table.
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
- `quantification_tool`: `"star"`, `"kallisto"`, or `"featurecounts"`.

## Tool Configuration

Thread and memory settings are in `config/tools.json`.

Important STAR options:

- `star.mem`: memory reserved by Snakemake.
- `star.threads`: threads passed to STAR.
- `star.read_files`: command used by STAR to read compressed FASTQ files, for
  example `gunzip -c`.
- `star.type`: STAR BAM output type, for example `BAM Unsorted`.
- `star.quant_mode`: STAR quantification mode used for transcriptome output.
- `star.gene_counts_strandedness`: count column used from
  `ReadsPerGene.out.tab`.
- `star.other`: extra STAR arguments.

Kallisto options:

- `kallisto.mem`: memory reserved by Snakemake.
- `kallisto.threads`: threads passed to Kallisto.

FeatureCounts options:

- `featurecounts.threads`: threads passed to FeatureCounts.
- `featurecounts.mem`: memory reserved by Snakemake.
- `featurecounts.strandedness`: FeatureCounts `-s` value. Use `0` for
  unstranded, `1` for stranded, or `2` for reversely stranded.
- `featurecounts.feature_type`: annotation feature type passed with `-t`, for
  example `exon` or `CDS`.
- `featurecounts.attribute`: annotation attribute passed with `-g`, for example
  `gene_id`.
- `featurecounts.extra`: additional FeatureCounts arguments.

QC and report options:

- `fastqc.mem` and `fastqc.threads`
- `fastq_screen.mem`, `fastq_screen.threads`, `fastq_screen.subset`, and
  `fastq_screen.aligner`
- `trim_galore.mem`
- `multiqc.mem`
- `gunc.mem` and `gunc.threads`

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

To run a specific quantification mode without editing `config/config.json`:

```bash
snakemake --use-conda --cores 4 --config quantification_tool=featurecounts
snakemake --use-conda --cores 4 --config quantification_tool=kallisto
snakemake --use-conda --cores 4 --config quantification_tool=star
```

On SLURM with Snakemake's executor support, adapt this pattern to your cluster:

```bash
snakemake --use-conda --slurm --default-resources slurm_account=YOUR_ACCOUNT --jobs 50
```

## Common Workflows

### Run Alignment-Based Gene Counts With STAR

Set `quantification_tool` to `star`, then run:

```bash
snakemake --use-conda --cores 4
```

Use this mode when you want STAR gene counts and STAR mapping metrics directly.

### Run Alignment-Based Gene Counts With FeatureCounts

Set `quantification_tool` to `featurecounts`, adjust the FeatureCounts settings
in `config/tools.json`, then run:

```bash
snakemake --use-conda --cores 4
```

Use this mode when you want explicit control over annotation feature type,
attribute, strandedness, and extra FeatureCounts options.

### Run Transcript Quantification With Kallisto

Set `quantification_tool` to `kallisto`, then run:

```bash
snakemake --use-conda --cores 4
```

Use this mode when transcript-level estimated counts and TPMs are desired.

## Current Scope And Limits

The workflow currently provides a basic DESeq2 analysis rather than a fully
general statistical modeling interface.

Current behavior:

- One metadata variable is used in the DESeq2 design.
- One configured reference level is used for the DESeq2 comparison.
- Count and metadata sample IDs must match exactly.
- Kallisto transcript counts are passed into the same DESeq2 entry point as the
  gene count matrices.
- Plots are generated for the configured DESeq2 comparison.

Common extensions that are not currently implemented include multi-factor
design formulas, explicit contrast tables, PCA plots, sample distance heatmaps,
MA plots, and automated replicate checks.

## Common Problems

- **Conda solver conflicts**: use the provided environment YAMLs with
  `conda-forge`, `bioconda`, and `defaults` channel priority. If the outer
  `snakemake` environment prints a `conda-libmamba-solver` warning, that is an
  issue with the active Snakemake environment rather than a workflow rule.
- **Metadata/count mismatch in DESeq2**: sample IDs in the metadata must match
  count-matrix columns exactly. The workflow stops with a clear error when they
  do not match.
- **Paired-end FASTQs not found**: check R1/R2 naming. The workflow recognizes
  common endings such as `_R1`, `_R2`, `_1.fastq`, and `_2.fastq`.
- **Wrong STAR count column**: set `star.gene_counts_strandedness` to
  `unstranded`, `forward`, or `reverse` to match the library strandedness.
- **Unexpected FeatureCounts assignment rates**: check
  `featurecounts.strandedness`, `featurecounts.feature_type`, and
  `featurecounts.attribute`.
- **Run location**: run the main workflow from the repository root, not from
  `test/` or `workflow/`.

## Development Checks

Use the Snakemake Conda environment from the repository root:

```bash
conda activate snakemake
snakemake --use-conda --cores 1 -n
snakemake --lint
```

The lint command can report style warnings for hardcoded output directories and
helper functions inside rule files. Treat DAG failures, environment solve
failures, and script errors as higher priority than style warnings.

## DAG

An example workflow DAG is available at:

![Workflow DAG](images/dag.svg)
