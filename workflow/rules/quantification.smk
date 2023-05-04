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
        reads=rules.trim_galore.output
    output:"results/{project}/quantification/kallisto/{subsample}/abundance.h5"
    log: "logs/{project}/quantification/kallisto/{subsample}.log"
    params:
        out_dir="results/{project}/quantification/kallisto/{subsample}/"
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        mkdir -p {params.out_dir}
        kallisto quant -i resources/Robinsoniella_peoriensis_annotation.feature_dna.fasta.transcriptome_index -o {params.out_dir} {input.reads}
        """

rule merge_kallisto:
    input: 
    output:
        "results/{project}/quantification/kallisto/transcript_tpms_all_samples.tsv",
        "results/{project}/quantification/kallisto/transcript_counts_all_samples.tsv"
    params:
        out_dir="results/{project}/quantification/kallisto/"
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        chmod 777 workflow/scripts/combine.sh
        workflow/scripts/combine.sh -i {params.out_dir} -o1 {output[0]}  -o2 {output[1]}
        """
