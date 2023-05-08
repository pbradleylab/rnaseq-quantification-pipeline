""" Add rules to this section that are related to alignment. This can be specific bam/sam processing steps 
or alternative aligners that may be run in leu of the default aligner hisat2.
"""
include: "trimming.smk"
include: "resources.smk"
configfile: "config/config.json"

# Compared to STAR, Hisat2 is more splice-junction aware aligner. Most bacterial RNA 
# transcripts do not undergo splicing; however, self-splicing introns may exist. Regardless,
# Hisat2 is a good choice for bacteria alignment even if no-splicing occurs since the absence
# will not affect the end results.
# Read more at the offical documentation for the tool: http:/daehwankimlab.github.io/hisat2/
rule hisat2:
    input: 
        rules.trim_galore.output,
        rules.hisat2_index.output
    output: "results/{project}/alignment/hisat2/{subsample}.bam"
    params:
        genome=config["genome"]
    log: "logs/{project}/alignment/hisat2/{subsample}.log"
    resources: mem=config["hisat2"]["mem"]
    conda: "../envs/alignment.yml"
    shell:
        """
        reads=$(echo {input[0]} | sed 's/ /,/g')
        hisat2 -U $reads -x {params.genome} | samtools sort -o {output}
        """

# Identifies duplicate reads to reduce over-representation.
rule mark_duplicates:
    input: 
        bam = rules.hisat2.output[0],
        bai = rules.hisat2.output[0]+".bai"
    output:
        bam="results/{project}/alignment/mark_duplicates/{subsample}.bam",
        metrics="results/{project}/alignment/mark_duplicates/{subsample}.metrics.txt"
    log: "logs/{project}/alignment/mark_duplicates/{subsample}.log"
    params:
        extra=config["picard"]["mark_duplicates"]["extra"],
        java_flags = lambda wildcards, input, output, threads, resources: config["java_flags_inject"].format(mem=resources.mem_mb, tmp=resources.tmpdir)
    resources:mem_mb=config["picard"]["mark_duplicates"]["mem"]
    conda: "../envs/alignment.yml"
    shell:"{params.java_flags} picard MarkDuplicates {params.extra} I={input.bam} O={output.bam} m={output.metrics} &> {log}"

# Calculate coverage per gene.
rule coverage_per_exon:
    input: 
        bam = rules.mark_duplicates.output[0]
    output: "{project}/reports/coverage_per_exon/{subsample}.coverage"
    params:
        subsample="{subsample}"
    log: "logs/{project}/reports/coverage_per_exon/{subsample}.log"
    conda: "../envs/quality_control.yml"
    shell:
        """
        bedtools bamtobed -i {input}| bedtools coverage -b - -a {params.subsample}.bed > {output}
        """ 