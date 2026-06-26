configfile: "config/config.json"
configfile: "config/tools.json"


def get_quant_method(wildcards):
    if config["quantification_tool"].lower() == "kallisto":
        return rules.merge_kallisto.output.counts.format(project=wildcards.project)
    if config["quantification_tool"].lower() == "star":
        return rules.merge_star_counts.output[0].format(project=wildcards.project)
    if config["quantification_tool"].lower() == "featurecounts":
        return rules.merge_featurecounts.output[0].format(project=wildcards.project)
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
        plot_svg="results/{project}/differential_expression/{project}.svg",
        plot_pdf="results/{project}/differential_expression/{project}.pdf",
        expression_boxplot="results/{project}/differential_expression/{project}_normalized_expression_boxplot.png",
        expression_boxplot_svg="results/{project}/differential_expression/{project}_normalized_expression_boxplot.svg",
        expression_boxplot_pdf="results/{project}/differential_expression/{project}_normalized_expression_boxplot.pdf",
        expression_density="results/{project}/differential_expression/{project}_normalized_expression_density.png",
        expression_density_svg="results/{project}/differential_expression/{project}_normalized_expression_density.svg",
        expression_density_pdf="results/{project}/differential_expression/{project}_normalized_expression_density.pdf",
        sample_distance_heatmap="results/{project}/differential_expression/{project}_sample_distance_heatmap.png",
        sample_distance_heatmap_svg="results/{project}/differential_expression/{project}_sample_distance_heatmap.svg",
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
        Rscript workflow/scripts/DESeq.R --counts_data {input.counts} --metadata_file {input.metadata} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.diffexp} --plot_path {output.plot} --plot_svg_path {output.plot_svg} --plot_pdf_path {output.plot_pdf} --expression_boxplot_path {output.expression_boxplot} --expression_boxplot_svg_path {output.expression_boxplot_svg} --expression_boxplot_pdf_path {output.expression_boxplot_pdf} --expression_density_path {output.expression_density} --expression_density_svg_path {output.expression_density_svg} --expression_density_pdf_path {output.expression_density_pdf} --sample_distance_heatmap_path {output.sample_distance_heatmap} --sample_distance_heatmap_svg_path {output.sample_distance_heatmap_svg} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """
