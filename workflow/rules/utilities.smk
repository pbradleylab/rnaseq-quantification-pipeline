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
        ln -s $PWD/{input.read1} {output.read1}
        ln -s $PWD/{input.read2} {output.read2}
        """

# Generates an bai index for a bam
rule samtools_index:
    input: "{BAM}.bam"
    output: "{BAM}.bam.bai"
    resources:mem_mb=1024
    conda: "../envs/alignment.yml"
    shell:"samtools index {input}"
