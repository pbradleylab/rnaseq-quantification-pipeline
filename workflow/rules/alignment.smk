""" Add rules to this section that are related to alignment. This can be specific bam/sam processing steps 
or alternative aligners that may be run in leu of the default aligner hisat2.
"""
include: "trimming.smk"
configfile: "config/config.json"

# Compared to STAR, Hisat2 is more splice-junction aware aligner. Most bacterial RNA 
# transcripts do not undergo splicing; however, self-splicing introns may exist. Regardless,
# Hisat2 is a good choice for bacteria alignment even if no-splicing occurs since the absence
# will not affect the end results.
# Read more at the offical documentation for the tool: http:/daehwankimlab.github.io/hisat2/
rule hisat2:
    input: rules.trim_galore.output
    output: "results/{project}/alignment/hisat2/{subsample}_se/{subsample}.sam"
    params:
        genome="resources/genomes/hg38"
    log: "logs/{project}/alignment/hisat2/{subsample}.log"
    resources: mem=config["hisat2"]["mem"]
    conda: "../envs/alignment.yml"
    shell:
        """
        reads=$(echo {input} | sed 's/ /,/g')
        hisat2 -U $reads -x {params.genome} -S {output} 2> {log}
        """