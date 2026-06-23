configfile: "config/config.json"
configfile: "config/tools.json"


def get_quant_method(wildcards):
    if config["quantification_tool"].lower() == "kallisto":
        return rules.merge_kallisto.output.counts.format(project=wildcards.project)
    if config["quantification_tool"].lower() == "star":
        return rules.merge_star_counts.output[0].format(project=wildcards.project)
    raise ValueError(
        f"Unsupported quantification_tool: {config['quantification_tool']}"
    )


rule deseq2:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        multiqc=rules.multiqc.output,
    output:
        diffexp="results/{project}/differential_expression/{project}.tsv",
        plot="results/{project}/differential_expression/{project}.png",
        plot_pdf="results/{project}/differential_expression/{project}.pdf",
    log:
        "logs/{project}/differential_expression/deseq2.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        reference_in_variable=config["deseq2"]["reference_in_variable"],
        variable_to_analyze=config["deseq2"]["variable_to_analyze"],
        log2fc_threshold=config["deseq2"].get("log2fc_threshold", 0.6),
        padj_threshold=config["deseq2"].get("padj_threshold", 0.05),
        label_top_n=config["deseq2"].get("label_top_n", 10),
    shell:
        """
        Rscript workflow/scripts/DESeq.R --counts_data {input.counts} --metadata_file {input.metadata} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.diffexp} --plot_path {output.plot} --plot_pdf_path {output.plot_pdf} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """
