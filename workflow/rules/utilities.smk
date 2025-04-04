""" Add rules to this section that are general utility functions or specific utilities.
"""
def get_r1_fastqs(wildcards):
    out = {} # Dictionary that will hold two reads with r1 in index 0 and r2 in index 1
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    r1=[x for x in reads if ("_R1" in x or ".R1" in x or ".r1" in x or "_r1" in x or "_1.fq" in x or "_1.fastq" in x)]
    return r1[0]

def get_r2_fastqs(wildcards):
    out = {} # Dictionary that will hold two reads with r1 in index 0 and r2 in index 1
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    r2=[x for x in reads if ("_R2" in x or ".R2" in x or ".r2" in x or "_r2" in x or "_2.fq" in x or "_2.fastq" in x)]
    return r2[0]

# Symlink reference genome
rule symlink_genome:
    input: config["genome"]["path"]
    output: "resources/symlink_genome/"+config["genome"]["name"]+".fa"
    shell:
        """
        ln -s {input} {output}
        """

# Symlink reference genome annotations in gff3 format
rule symlink_gff3:
    input: config["gff3"]["path"]
    output: "resources/symlink_gff3/"+config["genome_name"]+".gff3"
    shell:
        """
        ln -s {input} {output}
        """


# Create symlink that has the file renamed to specific ending
rule symlink_files:
    input: 
        read1=get_r1_fastqs,
        read2=get_r2_fastqs
    output: 
        read1="resources/{project}/raw/{subsample}_R1.fastq.gz",
        read2="resources/{project}/raw/{subsample}_R2.fastq.gz"
    params:
        outdir="resources/{project}/raw/"
    resources:mem_mb=1024
    shell:
        """
        ln -s {input.read1} {output.read1}
        ln -s {input.read2} {output.read2}
        """
