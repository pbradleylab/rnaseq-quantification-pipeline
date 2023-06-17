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

rule star_reads_per_gene:
    input: 
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.star_index.output
    output: "results/{project}/quantification/star/{subsample}/ReadsPerGene.out.tab"
    log: "logs/{project}/quantification/star_gene/{subsample}.log"
    params:
        outdir="results/{project}/quantification/star/{subsample}/",
        genome="resources/"+config["genome_name"]+"/"+config["genome_name"]+".transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=config["gff3"],
        bam_type=config["star"]["type"]
    threads:config["star"]["threads"]
    resources: mem=config["star"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.reads} \
                --readFilesCommand {params.read_files} --outSAMtype {params.bam_type} \
                --quantMode GeneCounts --outFileNamePrefix {params.outdir} %> {log}
        """

rule star_reads_per_transcript:
    input: 
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.star_index.output
    output: "results/{project}/quantification/star/{subsample}/Aligned.out.tab"
    log: "logs/{project}/quantification/star_transcript/{subsample}.log"
    params:
        outdir="results/{project}/quantification/star/{subsample}/",
        genome="resources/"+config["genome_name"]+"/"+config["genome_name"]+".transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=config["gff3"],
        bam_type=config["star"]["type"]
    threads:config["star"]["threads"]
    resources: mem=config["star"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.reads} \
                --readFilesCommand {params.read_files} --outSAMtype {params.bam_type} \
                --quantMode TranscriptomeSAM --outFileNamePrefix {params.outdir} %> {log}
        """