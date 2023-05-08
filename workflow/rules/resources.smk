""" Add rules to this section that are related to resource tranfsormtion for use in workflow, 
retrival, and mangement.
"""

def get_r1_fastqs(wildcards):
    out = {} # Dictionary that will hold two reads with r1 in index 0 and r2 in index 1
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    r1=[x for x in reads if ("_R1" in x or ".R1" in x or ".r1" in x or "_r1" in x or "_2.fq" in x or "_1.fastq" in x)]
    return r1[0]

def get_r2_fastqs(wildcards):
    out = {} # Dictionary that will hold two reads with r1 in index 0 and r2 in index 1
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    r2=[x for x in reads if ("_R2" in x or ".R2" in x or ".r2" in x or "_r2" in x or "_2.fq" in x or "_2.fastq" in x)]
    return r2[1]

# Create symlink that has the file renamed to specific ending
rule symlink_files:
    input: 
        read1=get_r1_fastqs,
        read2=get_r2_fastqs
    output: 
        read1="resources/{project}/raw/{subsample}_R1.fastq.gz",
        read2="resources/{project}/raw/{subsample}_R2.fastq.gz"
    params:
        out_dir="resources/{project}/raw/"
    resources:mem_mb=1024
    shell:
        """
        mkdir -p {params.out_dir}
        ln -s {input.read1} {output.read1}
        ln -s {input.read2} {output.read2}
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

rule hisat2_index:
    input: config["genome"]
    output: "resources/"+config["genome_name"]+".1.ht2"
    params:
        name=config["genome_name"]
    conda: "../envs/alignment.yml"
    shell:
        """
        hisat2-build {input} {params.name}
        mv {params.name}.*.ht2 resources/
        """     

# Generates an bai index for a bam
rule samtools_index:
    input: "{subsample}.bam"
    output: "{subsample}.bam.bai"
    resources:mem_mb=1024
    conda: "../envs/alignment.yml"
    shell:"samtools index {input}"