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

# Downlaod the default databases fro a gunc run on input reference
rule download_gunc_db:
    output: 
        output_dir=directory("resources/gunc/"),
        db="resources/gunc/gunc_db_progenomes2.1.dmnd"
    conda: "../envs/resources.yml"
    shell:
        """
        mkdir -p {output}
        gunc download_db {output}
        """

# Generate the transcriptome file.
rule gffread:
    input:
        gff3=config["gff3"],
        ref=config["genome"]
    output: "resources/"+config["genome_name"]+"/"+config["genome_name"]+".fa"
    conda: "../envs/resources.yml"
    shell:
        """
        gffread -w {output} -g {input.ref} --gtf {input.gff3}
        """

# Generate the Kallisto index
rule kallisto_index:
    input: rules.gffread.output
    output: "resources/"+config["genome_name"]+"/"+config["genome_name"]+".kallisto.index"
    params:
        name=config["genome_name"]
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        kallisto index -i {params.name}.kallisto.index {input} -h
        mv {params.name}.kallisto.index resources/{params.name}
        """

# Generate the STAR index
rule star_index:
    input: rules.gffread.output
    output: "resources/"+config["genome_name"]+"/"+config["genome_name"]+".kallisto.index"
    params:
        genome_dir=config["star"]["genome_dir"]
    resources: mem=config["star"]["mem"]
    threads: config["star"]["threads"]
    conda: "../envs/quantification.yml"
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.genome_dir} --genomeFastaFiles {input}
        """