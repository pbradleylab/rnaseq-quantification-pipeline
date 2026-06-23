configfile: "config/config.json"
configfile: "config/tools.json"

def get_quant_method(wildcards):
    if config["quantification_tool"].lower() == "kallisto":
        return rules.merge_kallisto.output.counts.format(project=wildcards.project)
    if config["quantification_tool"].lower() == "star":
        return rules.merge_star_counts.output[0].format(project=wildcards.project)
    raise ValueError(f"Unsupported quantification_tool: {config['quantification_tool']}")


rule deseq2:
    input:
        counts=get_quant_method,
        metadata=config["metadata"],
        multiqc=rules.multiqc.output
    output:
        diffexp="results/{project}/differential_expression/{project}.tsv",
        plot="results/{project}/differential_expression/{project}.png"
    params:
        reference_in_variable=config["deseq2"]["reference_in_variable"],
        variable_to_analyze=config["deseq2"]["variable_to_analyze"]       
    log: "logs/{project}/differential_expression/deseq2.log"
    conda: "../envs/DESeq2.yml"
    shell:
        """
        Rscript workflow/scripts/DESeq.R --counts_data {input.counts} --metadata_file {input.metadata} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.diffexp} --plot_path {output.plot} 2> {log}
        """
