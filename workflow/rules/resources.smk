""" Add rules to this section that are related to resource tranfsormtion for use in workflow, 
retrival, and mangement.
"""
# Downloads the default necessary resources for checking for contamination 
# with fastqscreen. Admittedly, not all of the genomes need/are to be used.
rule download_fastq_screen_genomes:
    output: directory("resources/FastQ_Screen_Genomes")
    params:
        folder_name=directory("FastQ_Screen_Genomes")
    log: "logs/resources/download_fastq_screen_genomes.log"
    conda: "../envs/resources.yml"
    shell: 
        """
        tmpdir=$(mktemp -d resources/.fastq_screen_genomes.XXXXXX)
        patched_fastq_screen="$tmpdir/fastq_screen"
        cp "$(command -v fastq_screen)" "$patched_fastq_screen"
        python3 workflow/scripts/patch_fastq_screen_get_genomes.py "$patched_fastq_screen"
        chmod +x "$patched_fastq_screen"
        (
            cd "$tmpdir"
            ./fastq_screen --get_genomes
        ) > {log} 2>&1
        rm -rf {output}
        mv "$tmpdir/{params.folder_name}" {output} >> {log} 2>&1
        rm -rf "$tmpdir"
        """

# Download given genome if needed
rule download_genome:
    output: "resources/download_genome/"
    params:
        url=config["genome"]["url"]
    log: "logs/resources/download_genome.log"
    conda: "../envs/resources.yml"
    shell:
        """
        wget -c --no-http-keep-alive {params.url} -O {output} > {log} 2>&1
        """

# Download given gff3 if needed
rule download_gff3:
    output: "resources/download_gff3/"
    params:
        url=config["gff3"]["url"]
    log: "logs/resources/download_gff3.log"
    conda: "../envs/resources.yml"
    shell:
        """
        wget -c --no-http-keep-alive {params.url} -O {output} > {log} 2>&1
        """

# Downlaod the default databases fro a gunc run on input reference
rule download_gunc_db:
    output: 
        output_dir=directory("resources/gunc/"),
        db="resources/gunc/gunc_db_progenomes2.1.dmnd"
    log: "logs/resources/download_gunc_db.log"
    conda: "../envs/resources.yml"
    shell:
        """
        mkdir -p {output.output_dir}
        gunc download_db {output.output_dir} > {log} 2>&1
        """

# Generate the transcriptome file.
rule gffread:
    input:
        gff3=lambda wildcards: rules.symlink_gff3.output if config["gff3"]["is_local"] else rules.download_gff3.output,
        ref=lambda wildcards: rules.symlink_genome.output if config["genome"]["is_local"] else rules.download_genome.output
    output: f"resources/{config['genome']['name']}/{config['genome']['name']}.fa"
    log: "logs/resources/gffread.log"
    conda: "../envs/resources.yml"
    shell:
        """
        gffread -w {output} -g {input.ref} --gtf {input.gff3} -F 2> {log}
        """
        
# Generate the Kallisto index
rule kallisto_index:
    input: rules.gffread.output
    output: f"resources/{config['genome']['name']}/{config['genome']['name']}.kallisto.index"
    params:
        name=config["genome"]["name"]
    resources: mem_mb=config["kallisto"]["mem"]
    log: "logs/resources/kallisto_index.log"
    conda: "../envs/quantification.yml"
    shell:
        """
        kallisto index -i {params.name}.kallisto.index {input} -h > {log} 2>&1
        mv {params.name}.kallisto.index resources/{params.name} >> {log} 2>&1
        """

rule convert_gff_to_gtf:
    input: lambda wildcards: rules.symlink_gff3.output if config["gff3"]["is_local"] else rules.download_gff3.output
    output: f"resources/{config['genome']['name']}/star.gtf"
    log: "logs/resources/convert_gff_to_gtf.log"
    conda: "../envs/resources.yml"
    shell:
        """
        gffread {input} -T -o {output} 2> {log}
        """

# Generate the STAR index
rule star_index:
    input: 
        ref=lambda wildcards: rules.symlink_genome.output if config["genome"]["is_local"] else rules.download_genome.output,
        gff=rules.convert_gff_to_gtf.output
    output: directory("resources/star/")
    resources: mem_mb=config["star"]["mem"]
    threads: config["star"]["threads"]
    log: "logs/resources/star_index.log"
    conda: "../envs/quantification.yml"
    shell:
        """
        mkdir -p {output}
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gff} > {log} 2>&1
        """
