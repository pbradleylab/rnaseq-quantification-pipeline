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

# Download given genome if needed
rule download_genome:
    output: "resources/download_genome/"
    params:
        url=config["genome"]["url"]
    conda: "../envs/setup.yml"
    shell:
        """
        wget -c --no-http-keep-alive {params.url} -O {output}
        """

# Download given gff3 if needed
rule download_gff3:
    output: "resources/download_gff3/"
    params:
        url=config["gff3"]["url"]
    conda: "../envs/setup.yml"
    shell:
        """
        wget -c --no-http-keep-alive {params.url} -O {output}
        """

# Downlaod the default databases fro a gunc run on input reference
rule download_gunc_db:
    output: 
        output_dir=directory("resources/gunc/"),
        db="resources/gunc/gunc_db_progenomes2.1.dmnd"
    conda: "../envs/resources.yml"
    shell:
        """
        mkdir -p {output.output_dir}
        gunc download_db {output.output_dir}
        """

# Generate the transcriptome file.
rule gffread:
    input:
        gff3=lambda wildcards: rules.symlink_gff3.output if config["gff3"]["is_local"] else rules.download_gff3.output,
        ref=lambda wildcards: rules.symlink_genome.output if config["genome"]["is_local"] else rules.download_genome.output
    output: "resources/"+config["genome"]["name"]+"/"+config["genome"]["name"]+".fa"
    conda: "../envs/resources.yml"
    shell:
        """
        gffread -w {output} -g {input.ref} --gtf {input.gff3} -F
        """
        
# Generate the Kallisto index
rule kallisto_index:
    input: rules.gffread.output
    output: "resources/"+config["genome"]["name"]+"/"+config["genome"]["name"]+".kallisto.index"
    params:
        name=config["genome"]["name"]
    resources: mem=config["kallisto"]["mem"]
    conda: "../envs/quantification.yml"
    shell:
        """
        kallisto index -i {params.name}.kallisto.index {input} -h
        mv {params.name}.kallisto.index resources/{params.name}
        """

rule convert_gff_to_gtf:
    input: lambda wildcards: rules.symlink_gff3.output if config["gff3"]["is_local"] else rules.download_gff3.output
    output: "resources/"+config["genome_name"]+"/star.gtf"
    conda: "../envs/resources.yml"
    shell:
        """
        gffread {input} -T -o {output}
        """

# Generate the STAR index
rule star_index:
    input: 
        ref=lambda wildcards: rules.symlink_genome.output if config["genome"]["is_local"] else rules.download_genome.output,
        gff=rules.convert_gff_to_gtf.output
    output: directory("resources/star/")
    resources: mem=config["star"]["mem"]
    threads: config["star"]["threads"]
    conda: "../envs/quantification.yml"
    shell:
        """
        mkdir -p {output}
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gff}
        """
