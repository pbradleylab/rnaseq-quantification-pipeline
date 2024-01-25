configfile: "config/config.json"
configfile: "config/tools.json"

def get_quant_method(wildcards):
    method = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        if config["quantification_tool"].lower() == "kallisto":
            method.append(rules.merge_kallisto.output[1].format(project=project))
        elif config["quantification_tool"].lower() == "star":
            method.append(rules.star_reads_per_gene.output[0].format(project=project, subsample=subsample))
            method.append(rules.star_reads_per_transcript.output[0].format(project=project, subsample=subsample))

    return list(set(method))


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
    conda: "../envs/DESeq2.yml"
    shell:
        """
        Rscript workflow/scripts/DESeq.R --counts_data {input.counts} --metadata_file {input.metadata} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output.diffexp} --plot_path {output.plot}
        """
