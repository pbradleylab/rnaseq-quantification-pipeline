""" Add rules to this section that are related to quantification.
"""
def get_h5(wildcards):
    out = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        seq_method = get_seq_method(subsample)
        if seq_method == "paired_end":
            out.append(rules.kallisto.output[0].format(project=project, subsample=subsample))
        elif seq_method == "single_end":
            out.append(rules.kallisto_single.output[0].format(project=project, subsample=subsample))
    return out

def get_star_gene_counts(wildcards):
    out = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        if project != wildcards.project:
            continue
        seq_method = get_seq_method(subsample)
        if seq_method == "paired_end":
            out.append(rules.star_reads_per_gene.output[0].format(project=project, subsample=subsample))
        elif seq_method == "single_end":
            out.append(rules.star_reads_per_gene_single.output[0].format(project=project, subsample=subsample))
    return out

# kallisto is a program for quantifying abundances of transcripts from bulk 
# and single-cell RNA-Seq data, or more generally of target sequences using 
# high-throughput sequencing reads.
# Official documentation can be found here: https:/pachterlab.github.io/kallisto/
rule kallisto:
    input: 
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.kallisto_index.output
    output: "results/{project}/quantification/kallisto/{subsample}/abundance.h5",
    log: "logs/{project}/quantification/kallisto/{subsample}.log"
    params:
        outdir="results/{project}/quantification/kallisto/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index"
    threads:config["kallisto"]["threads"]
    resources: mem_mb=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        kallisto quant -t {threads} -i {input.index} -o {params.outdir} {input.reads}
        """

rule kallisto_single:
    input:
        reads=rules.trim_galore_single.output,
        transcriptome=rules.gffread.output,
        index=rules.kallisto_index.output
    output:
        "results/{project}/quantification/kallisto_single/{subsample}/abundance.h5"
    log:
        "logs/{project}/quantification/kallisto_single/{subsample}.log"
    params:
        outdir="results/{project}/quantification/kallisto_single/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        frag_len=200,  # Set fragment length estimate
        frag_sd=20     # Set fragment length standard deviation
    threads: config["kallisto"]["threads"]
    resources:
        mem_mb=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        kallisto quant -t {threads} -i {input.index} -o {params.outdir} --single \
        -l {params.frag_len} -s {params.frag_sd} {input.reads} 2> {log}
        """
        
rule merge_kallisto:
    input: get_h5
    output:
        tpm="results/{project}/final/quantification/kallisto/transcript_tpms_all_samples.tsv",
        counts="results/{project}/final/quantification/kallisto/transcript_counts_all_samples.tsv"
    params:
        outdir="results/{project}/quantification/kallisto/",
        project="{project}"
    resources: mem_mb=config["kallisto"]["mem"]
    log: "logs/{project}/quantification/kallisto/merge.log"
    conda: "../envs/quantification.yml"
    shell:
        """
        bash workflow/scripts/combine_count.sh {params.project} {params.outdir} {output.tpm} 2> {log}
        bash workflow/scripts/combine_tpm.sh {params.project} {params.outdir} {output.counts} 2>> {log}
        """

rule star_reads_per_gene:
    input: 
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.star_index.output
    output: "results/{project}/quantification/star/star_reads_per_gene/{subsample}/ReadsPerGene.out.tab"
    log: "logs/{project}/quantification/star_gene/{subsample}.log"
    params:
        outdir="results/{project}/quantification/star/star_reads_per_gene/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=lambda wildcards: rules.symlink_gff3.output if config["gff3"]["is_local"] else rules.download_gff3.output,
        bam_type=config["star"]["type"],
        other=config["star"]["other"]
    threads:config["star"]["threads"]
    resources: mem_mb=config["star"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.reads} \
                --readFilesCommand {params.read_files} --outSAMtype {params.bam_type} \
                --quantMode GeneCounts --outFileNamePrefix {params.outdir} 2> {log} \
                {params.other}
        """

rule star_reads_per_gene_single:
    input:
        reads=rules.trim_galore_single.output,
        transcriptome=rules.gffread.output,
        index=rules.star_index.output
    output:
        "results/{project}/quantification/star/star_reads_per_gene_single/{subsample}/ReadsPerGene.out.tab"
    log:
        "logs/{project}/quantification/star_gene_single/{subsample}.log"
    params:
        outdir="results/{project}/quantification/star/star_reads_per_gene_single/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=lambda wildcards: rules.symlink_gff3.output if config["gff3"]["is_local"] else rules.download_gff3.output,
        bam_type=config["star"]["type"],
        other=config["star"]["other"]
    threads: config["star"]["threads"]
    resources:
        mem_mb=config["star"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.reads} \
             --readFilesCommand {params.read_files} --outSAMtype {params.bam_type} \
             --quantMode GeneCounts --outFileNamePrefix {params.outdir} \
             --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 2> {log} \
             {params.other}
        """

rule star_reads_per_transcript:
    input: 
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.star_index.output
    output: "results/{project}/quantification/star/star_reads_per_transcript/{subsample}/Aligned.toTranscriptome.out.bam"
    log: "logs/{project}/quantification/star_transcript/{subsample}.log"
    params:
        outdir="results/{project}/quantification/star/star_reads_per_transcript/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=lambda wildcards: rules.symlink_gff3.output if config["gff3"]["is_local"] else rules.download_gff3.output,
        bam_type=config["star"]["type"],
        other=config["star"]["other"]
    threads:config["star"]["threads"]
    resources: mem_mb=config["star"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.reads} \
                --readFilesCommand {params.read_files} --outSAMtype {params.bam_type} \
                --quantMode TranscriptomeSAM --outFileNamePrefix {params.outdir} 2> {log} \
                {params.other}
        """

rule star_reads_per_transcript_single:
    input:
        reads=rules.trim_galore_single.output,
        transcriptome=rules.gffread.output,
        index=rules.star_index.output
    output:
        "results/{project}/quantification/star/star_reads_per_transcript_single/{subsample}/Aligned.toTranscriptome.out.bam"
    log:
        "logs/{project}/quantification/star_transcript_single/{subsample}.log"
    params:
        outdir="results/{project}/quantification/star/star_reads_per_transcript_single/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=lambda wildcards: rules.symlink_gff3.output if config["gff3"]["is_local"] else rules.download_gff3.output,
        bam_type=config["star"]["type"],
        other=config["star"]["other"]
    threads: config["star"]["threads"]
    resources:
        mem_mb=config["star"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.reads} \
             --readFilesCommand {params.read_files} --outSAMtype {params.bam_type} \
             --quantMode TranscriptomeSAM --outFileNamePrefix {params.outdir} \
             --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 2> {log} \
             {params.other}
        """

rule merge_star_counts:
    input: get_star_gene_counts
    output: "results/{project}/final/quantification/star/gene_counts_all_samples.tsv"
    resources: mem_mb=config["star"]["mem"]
    log: "logs/{project}/quantification/star/merge_counts.log"
    conda: "../envs/quantification.yml"
    shell:
        """
        python workflow/scripts/combine_star_counts.py --output {output} {input} 2> {log}
        """
