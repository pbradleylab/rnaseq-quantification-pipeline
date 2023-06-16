""" Add rules to this section that are related to quantification.
"""
include: "metrics.smk"
include: "trimming.smk"

# kallisto is a program for quantifying abundances of transcripts from bulk 
# and single-cell RNA-Seq data, or more generally of target sequences using 
# high-throughput sequencing reads.
# Official documentation can be found here: https:/pachterlab.github.io/kallisto/
rule kallisto:
    input: 
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.kallisto_index.output
    output: "results/{project}/quantification/kallisto/{subsample}/abundance.h5"
    log: "logs/{project}/quantification/kallisto/{subsample}.log"
    params:
        outdir="results/{project}/quantification/kallisto/{subsample}/",
        genome="resources/"+config["genome_name"]+"/"+config["genome_name"]+".transcriptome_index"
    threads:config["kallisto"]["threads"]
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        kallisto quant -t {threads} -i {input.index} -o {params.outdir} {input.reads}
        """
        
rule merge_kallisto:
    input: rules.multiqc.output
    output:
        "results/{project}/final/quantification/kallisto/transcript_tpms_all_samples.tsv",
        "results/{project}/final/quantification/kallisto/transcript_counts_all_samples.tsv"
    params:
        outdir="results/{project}/quantification/kallisto/",
        project="{project}"
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        bash workflow/scripts/combine_count.sh {params.project} {params.outdir} {output[0]}
        bash workflow/scripts/combine_tpm.sh {params.project} {params.outdir} {output[1]}
        """

rule star:
    input: 
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.star_index.output
    output: "results/{project}/quantification/star/{subsample}/abundance.h5"
    params:
        outdir="results/{project}/quantification/star/{subsample}/",
        genome="resources/"+config["genome_name"]+"/"+config["genome_name"]+".transcriptome_index",
        gff=config["gff3"]
    threads:config["star"]["threads"]
    resources: mem=config["star"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        STAR --genomeDir {input.index} --outFileNamePrefix {params.outdir} --readFilesCommand gunzip -c --readFilesIn {input.reads}
        """

# rule star:
#     input: rules.trim_galore.output
#     output: "../results/{project}/alignment/star/{subsample}Aligned.toTranscriptome.out.bam"
#     params:
#         name="{subsample}",
#         outputDir="../results/{project}/alignment/star/",
#         readFilesCommand=config["star"]["read_files"],
#         genomeDir="resources/star/"
#     log: "logs/{project}/quantification/star/{subsample}.log"
#     threads: config["star"]["threads"]
#     resources: mem=config["star"]["mem"]
#     conda: "../envs/quantification.yml"
#     shell:
#         """
#         STAR --runThreadN {threads} --genomeDir {params.genomeDir} --readFilesIn {input.reads} --readFilesCommand {params.readFilesCommand} --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outFileNamePrefix {params.name} %> {log}
#         mkdir -p {params.outputDir}
#         mv $(pwd)/{params.name}* {params.outputDir}
#         """
