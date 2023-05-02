include: "trimming.smk"
configfile: "config/config.json"

# More Splice-junction aware aligner.
rule hisat2:
    input: rules.trim_galore.output
    output: "../results/{project}/alignment/hisat2/{subsample}_se/{subsample}.sam"
    params:
        genome="../resources/genomes/hg38"
    log: "../logs/{project}/alignment/hisat2/{subsample}.log"
    resources: mem=config["hisat2"]["mem"]
    conda: "../envs/alignment.yml"
    shell:
        """
        reads=$(echo {input} | sed 's/ /,/g')
        hisat2 -U $reads -x {params.genome} -S {output} 2> {log}
        """