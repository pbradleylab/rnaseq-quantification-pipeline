""" Add rules to this section that are related to quality control and post processing.
"""
from scripts.utils import *
configfile: "config/config.json"
include: "utilities.smk"
include: "resources.smk"

def get_fastq_input(wildcards):
    return get_subsample_attributes(wildcards.subsample, "reads", pep)

def get_multiqc_subsamples(wildcards):
    metricsLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)

        # Always run rules on the outside
        metricsLST.append(rules.fastqc.output[0].format(project=project, subsample=subsample))
        metricsLST.append(rules.fastq_screen.output[0].format(project=project, subsample=subsample))
        if config["quantification_tool"].lower() == "kallisto":
            metricsLST.append(rules.kallisto.output[0].format(project=project, subsample=subsample))
        elif config["quantification_tool"].lower() == "star":
            metricsLST.append(rules.star_reads_per_gene.output[0].format(project=project, subsample=subsample))
            metricsLST.append(rules.star_reads_per_transcript.output[0].format(project=project, subsample=subsample))
        if config["run_gunc"].lower() == "true":
            metricsLST.append(rules.gunc_plot.output[0].format(project=project, subsample=subsample))

    return metricsLST

# Run Gunc on the reference being used to make sure base assumptions are
# okay about reference being used before inferring anything else based on it
# This shoudl be run for both in house and genomes from sources like NCBI
rule gunc_run:
    input:
        ref=config["genome"],
        db=rules.download_gunc_db.output.db
    output: "results/{project}/reports/gunc_run/diamond_output/{project}.diamond.progenomes_2.1.out"
    log: "logs/{project}/reports/gunc_run/{project}_run.log"
    params:
        out_dir=directory("results/{project}/reports/gunc_run/")
    threads: config["gunc"]["threads"]
    resources: mem_mb = config["gunc"]["mem"]
    conda: "../envs/metrics.yml"
    shell: 
        """
        mkdir -p {params.out_dir}
        gunc run --threads {threads} --temp_dir /tmp/ -i {input.ref} --db_file {input.db} --out_dir {params.out_dir} 2> {log}
        """

# Plot Gunc output from run
rule gunc_plot:
    input: 
        ref=config["genome"],
        diamond=rules.gunc_run.output
    output:directory("results/{project}/reports/gunc_plot/")
    log: "logs/{project}/reports/gunc_plot/gunc_plot.log"
    conda: "../envs/metrics.yml"
    shell:
        """
        mkdir -p {output}
        gunc plot --diamond_file {input.diamond} --out_dir {output}
        """

# Calculates the contamination accross in the samples.
# The config must be supplied when calling this rule.
# Offical documentation: https:/www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
rule fastq_screen:
    input:
        refs=rules.download_fastq_screen_genomes.output,
        read1=rules.symlink_files.output[0],
        read2=rules.symlink_files.output[1]
    output:"results/{project}/reports/fastq_screen/{subsample}/{subsample}_R1_screen.html"
    params:
        fastq_screen_config="config/fastq_screen_config.conf",
        subset=config["fastq_screen"]["subset"],
        aligner=config["fastq_screen"]["aligner"],
        outdir="results/{project}/reports/fastq_screen/{subsample}/",
        subsample="{subsample}"
    log: "logs/{project}/reports/fastq_screen/{subsample}.log"
    threads: config["fastq_screen"]["threads"]
    resources: mem_mb = config["fastq_screen"]["mem"]
    conda: "../envs/metrics.yml"
    shell: 
        """
        fastq_screen --outdir {params.outdir} --force --aligner {params.aligner} --conf {params.fastq_screen_config} --subset {params.subset} --threads {threads} {input.read1} {input.read2} 2> {log}
        """

# FastQC is a quality control application for high throughput sequence data. 
# It reads in sequence data in a variety of formats and can either provide an interactive 
# application to review the results of several different QC checks, or create an HTML 
# based report which can be integrated into a pipeline.
#Offical Documentation: https:/www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc:
    input:
        read1=rules.symlink_files.output[0],
        read2=rules.symlink_files.output[1]
    output: "results/{project}/reports/fastqc/{subsample}_R1_fastqc.zip"
    log: "logs/{project}/reports/fastqc/{subsample}.log"
    params:
        outdir="results/{project}/reports/fastqc/",
        subsample="{subsample}"
    threads: config["fastqc"]["threads"]
    resources: mem_mb = config["fastqc"]["mem"]
    conda: "../envs/metrics.yml"
    shell: 
        """
        mkdir -p {params.outdir}
        fastqc {input.read1} {input.read2} --extract -t {threads} -o {params.outdir} 2> {log}
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
    conda: "../envs/metrics.yml"
    shell:
        """
        multiqc -f {params.search} --cl-config log_filesize_limit:2000000000 --outdir {params.output_dir} 2> {log}
        """
