include: "quantification.smk"
# Perform differential expression analysis at the gene-level with DESeq2
# Produce PCAs, heatmaps, volcano plots. Code taken from rna-seq-pop under
# MIT Liscensing: https://github.com/sanjaynagi/rna-seq-pop/tree/master
rule differential_gene_expression:
    input:
        counts=rules.kallisto.output,
        metadata=config["metadata"],
        genes2transcripts=config["genes2transcripts"]
    output:
        csvs=expand("results/{project}/differential_gene_expression/genediff/{comp}.csv", comp=config["contrasts"]),
        xlsx="results/{project}/differential_gene_expression/genediff/"+config["genome_name"]+"_diffexp.xlsx",
        pca="results/{project}/differential_gene_expression/counts/PCA.pdf",
        countStats="results/{project}/differential_gene_expression/counts/countStatistics.tsv",
        normCounts="results/{project}/differential_gene_expression/counts/normCounts.tsv",
    group:
        "diffexp"
    priority: 20
    conda:
        "../envs/differential_expression.yaml"
    log:
        "logs/{project}/differential_gene_expression/differential_gene_expression.log",
    params:
        DEcontrasts=config["contrasts"],
    script:
        "scripts/DESeq.R"