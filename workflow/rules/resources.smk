""" Add rules to this section that are related to resource tranfsormtion for use in workflow, 
retrival, and mangement.
"""

# Downloads the default necessary resources for checking for contamination 
# with fastqscreen. Admittedly, not all of the genomes need/are to be used.
rule download_fastq_screen_genomes:
    output: directory("resources/FastQ_Screen_Genomes")
    params:
        folder_name=directory("FastQ_Screen_Genomes")
    conda: "../envs/resources.yml"
    shell: 
        """
        fastq_screen --get_genomes
        mv {params.folder_name} {output}
        """

# Generate the Kallisto index
rule kallisto_index:
    output: "resources/"+config["genome_name"]+".index"
    params:
        ref=config["genome"],
        name=config["genome_name"]
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        kallisto index -i {params.name}.index {params.ref}
        mv {params.name}.index resources/
        """