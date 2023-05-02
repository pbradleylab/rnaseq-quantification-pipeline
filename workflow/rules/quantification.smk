""" Add rules to this section that are related to quantification.
"""
include: "alignment.smk"
include: "resources.smk"

# kallisto is a program for quantifying abundances of transcripts from bulk 
# and single-cell RNA-Seq data, or more generally of target sequences using 
# high-throughput sequencing reads.
# Official documentation can be found here: https:/pachterlab.github.io/kallisto/
rule kallisto:
    input: 
        reads=rules.trim_galore.output,
        index=rules.kallisto_index.output
    output:"results/{project}/quantification/kallisto/{subsample}/abundance.h5"
    log: "logs/{project}/quantification/kallisto/{subsample}.log"
    params:
        out_dir="results/{project}/quantification/kallisto/{subsample}/"
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        mkdir -p {params.out_dir}
        kallisto quant -i {input.index} -o {params.out_dir} {input.reads}
        """
