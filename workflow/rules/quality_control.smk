""" Add rules to this section that are related to quality control and post processing.
"""
from scripts.utils import *
configfile: "config/config.json"
include: "differential_expression.smk"

def get_fastq_input(wildcards):
    return get_subsample_attributes(wildcards.subsample, "reads", pep)

def get_multiqc_subsamples(wildcards):
    metricsLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        
        # Always run rules on the outside
        metricsLST.append(rules.fastqc.output[0].format(project=project, subsample=subsample))
        metricsLST.append(rules.fastq_screen.output[0].format(project=project, subsample=subsample))
        metricsLST.append(rules.kallisto.output[0].format(project=project, subsample=subsample))
        #metricsLST.append(rules.differential_gene_expression.output[0].format(project=project, subsample=subsample))

    return metricsLST

# Calculates the contamination accross in the samples.
# The config must be supplied when calling this rule.
# Offical documentation: https:/www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
rule fastq_screen:
    input: 
        unpack(get_fastq_input),
        rules.download_fastq_screen_genomes.output
    output:"results/{project}/reports/fastq_screen/{subsample}/FastQ_Screen_Genomes_screen.txt"
    params:
        fastq_screen_config="config/fastq_screen_config.conf",
        subset=config["fastq_screen"]["subset"],
        aligner=config["fastq_screen"]["aligner"],
        out_dir="results/{project}/reports/fastq_screen/{subsample}/",
        subsample="{subsample}"
    log: "logs/{project}/reports/fastq_screen/{subsample}.log"
    threads: config["fastq_screen"]["threads"]
    resources: mem_mb = config["fastq_screen"]["mem"]
    conda: "../envs/quality_control.yml"
    shell: 
        """
        fastq_screen --outdir {params.out_dir} --force --aligner {params.aligner} --conf {params.fastq_screen_config} --subset {params.subset} --threads {threads} {input} 2> {log}
        """

# FastQC is a quality control application for high throughput sequence data. 
# It reads in sequence data in a variety of formats and can either provide an interactive 
# application to review the results of several different QC checks, or create an HTML 
# based report which can be integrated into a pipeline.
#Offical Documentation: https:/www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc:
    input: unpack(get_fastq_input)
    output: "results/{project}/reports/fastqc/{subsample}_R1_001_fastqc.zip"
    log: "logs/{project}/reports/fastqc/{subsample}.log"
    params:
        out_dir="results/{project}/reports/fastqc/",
        subsample="{subsample}"
    threads: config["fastqc"]["threads"]
    resources: mem_mb = config["fastqc"]["mem"]
    conda: "../envs/quality_control.yml"
    shell: 
        """
        mkdir -p {params.out_dir}
        fastqc {input} --extract -t {threads} -o {params.out_dir} 2> {log}
        """

rule samtools_depth:
    input: rules.mark_duplicates.output.bam
    output: "results/{project}/reports/samtools_depth/{subsample}.coverage"
    params:
        out_dir="results/{project}/reports/samtools_depth/",
        subsample="{subsample}"
    log: "logs/{project}/samtools_depth/{subsample}.log"
    threads: config["samtools_depth"]["threads"]
    resources: mem_mb = config["samtools_depth"]["mem"]
    conda: "../envs/quality_control.yml"
    shell: 
        """
        samtools depth {input} > {output} 2> {log}
        """

# Generates html reports of analysis / QC done using MultiQC
# as shown here: https:/multiqc.info/
# Files that are too large may be skipped and this rule may
# fail even if the fastqc is successfully run. To adjust this
# the config file for fastqc may be adjusted.
rule multiqc:
    input:get_multiqc_subsamples
    output: 
        html="results/{project}/final/multiqc/multiqc_report.html"
    params:
        search=directory("results"),
        output_dir=directory("results/{project}/final/multiqc/")
    log: "logs/{project}/reports/multiqc/multiqc_report.log"
    resources: mem_mb = config["multiqc"]["mem"]
    conda: "../envs/quality_control.yml"
    shell:
        """
        multiqc -f {params.search} --cl-config log_filesize_limit:2000000000 --outdir {params.output_dir} 2> {log}
        """