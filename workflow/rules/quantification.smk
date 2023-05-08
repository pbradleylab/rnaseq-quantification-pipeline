""" Add rules to this section that are related to quantification.
"""
include: "alignment.smk"
include: "quality_control.smk"

# kallisto is a program for quantifying abundances of transcripts from bulk 
# and single-cell RNA-Seq data, or more generally of target sequences using 
# high-throughput sequencing reads.
# Official documentation can be found here: https:/pachterlab.github.io/kallisto/
rule kallisto:
    input: 
        reads=rules.trim_galore.output,
        transcriptome=rules.gffread.output,
        index=rules.kallisto_index.output
    output:"results/{project}/quantification/kallisto/{subsample}/abundance.h5"
    log: "logs/{project}/quantification/kallisto/{subsample}.log"
    params:
        out_dir="results/{project}/quantification/kallisto/{subsample}/",
        genome="resources/"+config["genome_name"]+".transcriptome_index"
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        mkdir -p {params.out_dir}
        echo {input.reads}
        kallisto quant -i {params.genome} -o {params.out_dir} {input.reads}
        """

rule merge_kallisto:
    input: rules.multiqc.output
    output:
        "results/{project}/final/quantification/kallisto/transcript_tpms_all_samples.tsv",
        "results/{project}/final/quantification/kallisto/transcript_counts_all_samples.tsv"
    params:
        out_dir="results/{project}/quantification/kallisto/"
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        chmod 777 workflow/scripts/combine.sh
        echo workflow/scripts/combine.sh -i {params.out_dir} -t {output[0]} -c {output[1]}
        bash workflow/scripts/combine.sh -i {params.out_dir} -t {output[0]} -c {output[1]}
        """
