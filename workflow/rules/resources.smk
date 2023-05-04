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

# Generate the transcriptome file.
rule gffread:
    input:
        gff3=config["gff3"],
        ref=config["genome"]
    output: "resources/transcriptome.fa"
    conda: "../envs/resources.yml"
    shell:
        """
        gffread -w {output} -g {input.ref} --gtf {input.gff3}
        """

# Generate the Kallisto index
rule kallisto_index:
    input: rules.gffread.output
    output: "resources/"+config["genome_name"]+".transcriptome_index"
    params:
        ref=config["genome"],
        name=config["genome_name"]
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        kallisto index -i {params.name}.transcriptome_index {input}
        mv {params.name}.transcriptome_index resources/
        """