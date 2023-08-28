configfile: "config/config.json"
configfile: "config/tools.json"

rule deseq2:
    input: 
        counts=rules.merge_kallisto.output.counts,
        metadata=config["metadata"]
    output: "results/{project}/differential_expression/{project}.tsv"
    params:
        reference_in_variable=config["deseq2"]["reference_in_variable"],
        variable_to_analyze=config["deseq2"]["variable_to_analyze"]
    conda: "../envs/DESeq.yml"
    shell:
        """
        Rscript workflow/scripts/DESeq.R --counts_data {input.counts} --metadata_file {input.metadata} --variable_to_analyze {params.variable_to_analyze} --reference_in_variable {params.reference_in_variable} --output_file {output}
        """