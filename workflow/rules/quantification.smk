"""Add rules to this section that are related to quantification."""


# kallisto is a program for quantifying abundances of transcripts from bulk
# and single-cell RNA-Seq data, or more generally of target sequences using
# high-throughput sequencing reads.
# Official documentation can be found here: https:/pachterlab.github.io/kallisto/
rule kallisto:
    input:
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.kallisto_index.output,
    output:
        "results/{project}/quantification/kallisto/{subsample}/abundance.h5",
    log:
        "logs/{project}/quantification/kallisto/{subsample}.log",
    conda:
        "../envs/quantification.yml"
    threads: config["kallisto"]["threads"]
    resources:
        mem_mb=config["kallisto"]["mem"],
    params:
        outdir="results/{project}/quantification/kallisto/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
    shell:
        """
        kallisto quant -t {threads} -i {input.index} -o {params.outdir} {input.reads}
        """


rule kallisto_single:
    input:
        reads=rules.trim_galore_single.output,
        transcriptome=rules.gffread.output,
        index=rules.kallisto_index.output,
    output:
        "results/{project}/quantification/kallisto_single/{subsample}/abundance.h5",
    log:
        "logs/{project}/quantification/kallisto_single/{subsample}.log",
    conda:
        "../envs/quantification.yml"
    threads: config["kallisto"]["threads"]
    resources:
        mem_mb=config["kallisto"]["mem"],
    params:
        outdir="results/{project}/quantification/kallisto_single/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        frag_len=200,  # Set fragment length estimate
        frag_sd=20,  # Set fragment length standard deviation
    shell:
        """
        kallisto quant -t {threads} -i {input.index} -o {params.outdir} --single \
        -l {params.frag_len} -s {params.frag_sd} {input.reads} 2> {log}
        """


rule merge_kallisto:
    input:
        lambda wildcards: get_kallisto_h5(wildcards, pep, rules),
    output:
        tpm="results/{project}/final/quantification/kallisto/transcript_tpms_all_samples.tsv",
        counts="results/{project}/final/quantification/kallisto/transcript_counts_all_samples.tsv",
    log:
        "logs/{project}/quantification/kallisto/merge.log",
    conda:
        "../envs/quantification.yml"
    resources:
        mem_mb=config["kallisto"]["mem"],
    params:
        outdirs=lambda wildcards: [
            f"results/{wildcards.project}/quantification/kallisto/",
            f"results/{wildcards.project}/quantification/kallisto_single/",
        ],
        project="{project}",
    shell:
        """
        python workflow/scripts/combine_kallisto.py --counts {output.counts} --tpm {output.tpm} {params.outdirs} 2> {log}
        """


rule star_reads_per_gene:
    input:
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.star_index.output,
    output:
        counts="results/{project}/quantification/star/star_reads_per_gene/{subsample}/ReadsPerGene.out.tab",
        log_final="results/{project}/quantification/star/star_reads_per_gene/{subsample}/Log.final.out",
    log:
        "logs/{project}/quantification/star_gene/{subsample}.log",
    conda:
        "../envs/quantification.yml"
    threads: config["star"]["threads"]
    resources:
        mem_mb=config["star"]["mem"],
    params:
        outdir="results/{project}/quantification/star/star_reads_per_gene/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=lambda wildcards: (
            rules.symlink_gff3.output
            if config["gff3"]["is_local"]
            else rules.download_gff3.output
        ),
        bam_type=config["star"]["type"],
        other=config["star"]["other"],
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
        index=rules.star_index.output,
    output:
        counts="results/{project}/quantification/star/star_reads_per_gene_single/{subsample}/ReadsPerGene.out.tab",
        log_final="results/{project}/quantification/star/star_reads_per_gene_single/{subsample}/Log.final.out",
    log:
        "logs/{project}/quantification/star_gene_single/{subsample}.log",
    conda:
        "../envs/quantification.yml"
    threads: config["star"]["threads"]
    resources:
        mem_mb=config["star"]["mem"],
    params:
        outdir="results/{project}/quantification/star/star_reads_per_gene_single/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=lambda wildcards: (
            rules.symlink_gff3.output
            if config["gff3"]["is_local"]
            else rules.download_gff3.output
        ),
        bam_type=config["star"]["type"],
        other=config["star"]["other"],
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
        index=rules.star_index.output,
    output:
        bam="results/{project}/quantification/star/star_reads_per_transcript/{subsample}/Aligned.toTranscriptome.out.bam",
        log_final="results/{project}/quantification/star/star_reads_per_transcript/{subsample}/Log.final.out",
    log:
        "logs/{project}/quantification/star_transcript/{subsample}.log",
    conda:
        "../envs/quantification.yml"
    threads: config["star"]["threads"]
    resources:
        mem_mb=config["star"]["mem"],
    params:
        outdir="results/{project}/quantification/star/star_reads_per_transcript/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=lambda wildcards: (
            rules.symlink_gff3.output
            if config["gff3"]["is_local"]
            else rules.download_gff3.output
        ),
        bam_type=config["star"]["type"],
        other=config["star"]["other"],
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
        index=rules.star_index.output,
    output:
        bam="results/{project}/quantification/star/star_reads_per_transcript_single/{subsample}/Aligned.toTranscriptome.out.bam",
        log_final="results/{project}/quantification/star/star_reads_per_transcript_single/{subsample}/Log.final.out",
    log:
        "logs/{project}/quantification/star_transcript_single/{subsample}.log",
    conda:
        "../envs/quantification.yml"
    threads: config["star"]["threads"]
    resources:
        mem_mb=config["star"]["mem"],
    params:
        outdir="results/{project}/quantification/star/star_reads_per_transcript_single/{subsample}/",
        genome=f"resources/{config['genome']['name']}/{config['genome']['name']}.transcriptome_index",
        read_files=config["star"]["read_files"],
        gff=lambda wildcards: (
            rules.symlink_gff3.output
            if config["gff3"]["is_local"]
            else rules.download_gff3.output
        ),
        bam_type=config["star"]["type"],
        other=config["star"]["other"],
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.reads} \
             --readFilesCommand {params.read_files} --outSAMtype {params.bam_type} \
             --quantMode TranscriptomeSAM --outFileNamePrefix {params.outdir} \
             --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 2> {log} \
             {params.other}
        """


rule merge_star_counts:
    input:
        lambda wildcards: get_star_gene_counts(wildcards, pep, rules),
    output:
        "results/{project}/final/quantification/star/gene_counts_all_samples.tsv",
    log:
        "logs/{project}/quantification/star/merge_counts.log",
    conda:
        "../envs/quantification.yml"
    resources:
        mem_mb=config["star"]["mem"],
    params:
        strandedness=config["star"].get("gene_counts_strandedness", "unstranded"),
    shell:
        """
        python workflow/scripts/combine_star_counts.py --output {output} --strandedness {params.strandedness} {input} 2> {log}
        """
