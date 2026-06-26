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


def get_annotation(wildcards):
    if config["gff3"]["is_local"]:
        return rules.symlink_gff3.output[0]
    return rules.download_gff3.output[0]


def deseq_common_params(config):
    return {
        "design_formula": config["deseq2"].get(
            "design_formula", f"~ {config['deseq2']['variable_to_analyze']}"
        ),
        "reference_in_variable": config["deseq2"]["reference_in_variable"],
        "variable_to_analyze": config["deseq2"]["variable_to_analyze"],
        "log2fc_threshold": config["deseq2"].get("log2fc_threshold", 0.6),
        "padj_threshold": config["deseq2"].get("padj_threshold", 0.05),
        "label_top_n": config["deseq2"].get("label_top_n", 10),
        "transform_method": config["deseq2"].get("transform_method", "vst"),
    }


rule count_annotation_overlap:
    input:
        counts=get_quant_method,
        annotation=get_annotation,
    output:
        "results/{project}/reports/count_annotation_overlap/{project}_count_annotation_overlap.tsv",
    log:
        "logs/{project}/reports/count_annotation_overlap/count_annotation_overlap.log",
    params:
        feature_type=config["featurecounts"].get("feature_type", "gene"),
        attribute=config["featurecounts"].get("attribute", "ID"),
    shell:
        """
        python3 workflow/scripts/check_count_annotation_overlap.py --counts {input.counts} --annotation {input.annotation} --feature-type {params.feature_type} --attribute {params.attribute} --output {output} 2> {log}
        """


rule gene_biotype_count_summary:
    input:
        counts=get_quant_method,
        annotation=get_annotation,
    output:
        "results/{project}/reports/gene_biotype_count_summary/{project}_gene_biotype_count_summary.tsv",
    log:
        "logs/{project}/reports/gene_biotype_count_summary/gene_biotype_count_summary.log",
    params:
        feature_type=config["featurecounts"].get("feature_type", "gene"),
        attribute=config["featurecounts"].get("attribute", "ID"),
    shell:
        """
        python3 workflow/scripts/summarize_gene_biotype_counts.py --counts {input.counts} --annotation {input.annotation} --feature-type {params.feature_type} --id-attribute {params.attribute} --output {output} 2> {log}
        """


rule deseq2:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        annotation_overlap=rules.count_annotation_overlap.output,
    output:
        diffexp="results/{project}/differential_expression/deseq2/{project}.tsv",
    log:
        "logs/{project}/deseq2/results.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode results --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.diffexp} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_significant_gene_tables:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        all="results/{project}/differential_expression/significant_gene_tables/{project}_significant_genes.tsv",
        up="results/{project}/differential_expression/significant_gene_tables/{project}_significant_upregulated_genes.tsv",
        down="results/{project}/differential_expression/significant_gene_tables/{project}_significant_downregulated_genes.tsv",
    log:
        "logs/{project}/differential_expression/significant_gene_tables/significant_gene_tables.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode significant_genes --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.all} --up_file {output.up} --down_file {output.down} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_volcano_plot:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        png="results/{project}/differential_expression/volcano_plot/{project}_volcano_plot.png",
        svg="results/{project}/differential_expression/volcano_plot/{project}_volcano_plot.svg",
    log:
        "logs/{project}/differential_expression/volcano_plot/volcano_plot.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode volcano --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.png} --svg_file {output.svg} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_normalized_counts:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        "results/{project}/differential_expression/normalized_counts/{project}_normalized_counts.tsv",
    log:
        "logs/{project}/differential_expression/normalized_counts/normalized_counts.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode normalized_counts --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output} --transform_method {params.transform_method} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_transformed_counts:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        "results/{project}/differential_expression/transformed_counts/{project}_transformed_counts.tsv",
    log:
        "logs/{project}/differential_expression/transformed_counts/transformed_counts.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode transformed_counts --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output} --transform_method {params.transform_method} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_cooks_reports:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        gene="results/{project}/differential_expression/cooks_reports/{project}_cooks_gene_report.tsv",
        sample="results/{project}/differential_expression/cooks_reports/{project}_cooks_sample_report.tsv",
    log:
        "logs/{project}/differential_expression/cooks_reports/cooks_reports.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode cooks_report --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.gene} --sample_report_file {output.sample} --transform_method {params.transform_method} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_ma_plot:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        png="results/{project}/differential_expression/ma_plot/{project}_ma_plot.png",
        svg="results/{project}/differential_expression/ma_plot/{project}_ma_plot.svg",
    log:
        "logs/{project}/differential_expression/ma_plot/ma_plot.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode ma --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.png} --svg_file {output.svg} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_expression_boxplot:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        png="results/{project}/differential_expression/normalized_expression_boxplot/{project}_normalized_expression_boxplot.png",
        svg="results/{project}/differential_expression/normalized_expression_boxplot/{project}_normalized_expression_boxplot.svg",
        pdf="results/{project}/differential_expression/normalized_expression_boxplot/{project}_normalized_expression_boxplot.pdf",
    log:
        "logs/{project}/differential_expression/normalized_expression_boxplot/normalized_expression_boxplot.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode expression_boxplot --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.png} --svg_file {output.svg} --pdf_file {output.pdf} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_expression_density:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        png="results/{project}/differential_expression/normalized_expression_density/{project}_normalized_expression_density.png",
        svg="results/{project}/differential_expression/normalized_expression_density/{project}_normalized_expression_density.svg",
        pdf="results/{project}/differential_expression/normalized_expression_density/{project}_normalized_expression_density.pdf",
    log:
        "logs/{project}/differential_expression/normalized_expression_density/normalized_expression_density.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode expression_density --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.png} --svg_file {output.svg} --pdf_file {output.pdf} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_sample_distance_heatmap:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        png="results/{project}/differential_expression/sample_distance_heatmap/{project}_sample_distance_heatmap.png",
        svg="results/{project}/differential_expression/sample_distance_heatmap/{project}_sample_distance_heatmap.svg",
    log:
        "logs/{project}/differential_expression/sample_distance_heatmap/sample_distance_heatmap.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode sample_distance_heatmap --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.png} --svg_file {output.svg} --transform_method {params.transform_method} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_pca:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        png="results/{project}/differential_expression/pca/{project}_pca.png",
        svg="results/{project}/differential_expression/pca/{project}_pca.svg",
    log:
        "logs/{project}/differential_expression/pca/pca.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode pca --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.png} --svg_file {output.svg} --transform_method {params.transform_method} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """


rule deseq2_library_size_factors:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        diffexp=rules.deseq2.output.diffexp,
    output:
        png="results/{project}/differential_expression/library_sizes_size_factors/{project}_library_sizes_size_factors.png",
        svg="results/{project}/differential_expression/library_sizes_size_factors/{project}_library_sizes_size_factors.svg",
    log:
        "logs/{project}/differential_expression/library_sizes_size_factors/library_sizes_size_factors.log",
    conda:
        "../envs/DESeq2.yml"
    params:
        **deseq_common_params(config)
    shell:
        """
        Rscript workflow/scripts/DESeq.R --mode library_size_factors --counts_data {input.counts} --metadata_file {input.metadata} --design_formula {params.design_formula:q} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.png} --svg_file {output.svg} --log2fc_threshold {params.log2fc_threshold} --padj_threshold {params.padj_threshold} --label_top_n {params.label_top_n} 2> {log}
        """
