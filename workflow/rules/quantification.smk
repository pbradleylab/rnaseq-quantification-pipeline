include: "alignment.smk"

rule kallisto:
    input: rules.trim_galore.output
    output:"../results/{project}/quantification/salmon/{subsample}.sf"
    log: "../logs/{project}/quantification/salmon/{subsample}.log"
    params:
        ref=config["genome"],
        trans="../results/{project}/quantification/salmon/{subsample}"
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        salmon quant -t {params.ref} -l A -a {input} -o {params.trans} --gcBias
        mv {params.trans}/quant.sf {output}
        rm -r {params.trans}
        """