""" Add rules to this section that are related to trimming.
"""
# Trim Galore is a wrapper around Cutadapt and FastQC 
# to consistently apply adapter and quality trimming to 
# FastQ files, with extra functionality for RRBS data.
# Official documentation can be found here: https:/github.com/FelixKrueger/TrimGalore
rule trim_galore:
    input: 
        read1=rules.symlink_files.output[0],
        read2=rules.symlink_files.output[1]
    output:
        r1="results/{project}/trimming/trim_galore/{subsample}_R1_val_1.fq.gz",
        r2="results/{project}/trimming/trim_galore/{subsample}_R2_val_2.fq.gz"
    params:
        outdir="results/{project}/trimming/trim_galore/" 
    log: "logs/{project}/trimming/trim_galore/{subsample}.log"
    resources:mem_mb=config["trim_galore"]["mem"]
    conda: "../envs/trimming.yml"
    shell: "trim_galore --paired {input.read1} {input.read2} -o {params.outdir} 2> {log}"

rule trim_galore_single:
    input:
        read1=rules.symlink_files_single.output
    output: "results/{project}/trimming/trim_galore/{subsample}_trimmed.fq.gz"
    params:
        outdir="results/{project}/trimming/trim_galore/"
    log: "logs/{project}/trimming/trim_galore/{subsample}.log"
    resources: mem_mb=config["trim_galore"]["mem"]
    conda: "../envs/trimming.yml"
    shell: "trim_galore {input.read1} -o {params.outdir} 2> {log}"


rule sortmerna:
    input:
        read1=rules.trim_galore.output.r1,
        read2=rules.trim_galore.output.r2,
    output:
        r1="results/{project}/trimming/sortmerna_paired/{subsample}_R1.non_rRNA.fastq.gz",
        r2="results/{project}/trimming/sortmerna_paired/{subsample}_R2.non_rRNA.fastq.gz",
        report="results/{project}/reports/sortmerna_paired/{subsample}.sortmerna.log",
    log:
        "logs/{project}/trimming/sortmerna/{subsample}.log"
    conda:
        "../envs/trimming.yml"
    threads: config["sortmerna"]["threads"]
    resources:
        mem_mb=config["sortmerna"]["mem"],
    params:
        refs=lambda wildcards: " ".join(
            f"--ref {ref}" for ref in as_list(config["sortmerna"]["refs"])
        ),
        workdir="results/{project}/trimming/sortmerna_paired/{subsample}.work/",
        outdir="results/{project}/trimming/sortmerna_paired/",
        reportdir="results/{project}/reports/sortmerna_paired/",
        extra=config["sortmerna"].get("extra", ""),
    shell:
        """
        rm -rf {params.workdir}
        mkdir -p {params.workdir} {params.outdir} {params.reportdir}
        sortmerna {params.refs} --reads {input.read1} --reads {input.read2} \
            --threads {threads} --workdir {params.workdir} --aligned rRNA_reads \
            --fastx --other non_rRNA_reads --paired_in --out2 {params.extra} \
            > {log} 2>&1
        mv {params.workdir}/non_rRNA_reads_fwd.f*q.gz {output.r1}
        mv {params.workdir}/non_rRNA_reads_rev.f*q.gz {output.r2}
        mv {params.workdir}/rRNA_reads.log {output.report}
        """


rule sortmerna_single:
    input:
        read1=rules.trim_galore_single.output,
    output:
        reads="results/{project}/trimming/sortmerna_single/{subsample}.non_rRNA.fastq.gz",
        report="results/{project}/reports/sortmerna_single/{subsample}.sortmerna.log",
    log:
        "logs/{project}/trimming/sortmerna/{subsample}.log"
    conda:
        "../envs/trimming.yml"
    threads: config["sortmerna"]["threads"]
    resources:
        mem_mb=config["sortmerna"]["mem"],
    params:
        refs=lambda wildcards: " ".join(
            f"--ref {ref}" for ref in as_list(config["sortmerna"]["refs"])
        ),
        workdir="results/{project}/trimming/sortmerna_single/{subsample}.work/",
        outdir="results/{project}/trimming/sortmerna_single/",
        reportdir="results/{project}/reports/sortmerna_single/",
        extra=config["sortmerna"].get("extra", ""),
    shell:
        """
        rm -rf {params.workdir}
        mkdir -p {params.workdir} {params.outdir} {params.reportdir}
        sortmerna {params.refs} --reads {input.read1} --threads {threads} \
            --workdir {params.workdir} --aligned rRNA_reads --fastx \
            --other non_rRNA_reads {params.extra} > {log} 2>&1
        mv {params.workdir}/non_rRNA_reads.f*q.gz {output.reads}
        mv {params.workdir}/rRNA_reads.log {output.report}
        """
