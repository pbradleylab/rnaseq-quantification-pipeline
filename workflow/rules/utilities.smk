""" Add rules to this section that are general utility functions or specific utilities.
"""
def get_fastqs(wildcards):
    out = {"r1": [], "r2": [], "single": []}
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    if isinstance(reads, str):
        reads = [reads]

    seq_method = get_subsample_attributes(wildcards.subsample, "seq_method", pep)
    if isinstance(seq_method, list):
        seq_method = seq_method[0]

    for read in reads:
        if seq_method == "paired_end":
            if any(x in read for x in ["_R1", ".R1", ".r1", "_r1", "_1.fq", "_1.fastq"]):
                out["r1"].append(read)
            elif any(x in read for x in ["_R2", ".R2", ".r2", "_r2", "_2.fq", "_2.fastq"]):
                out["r2"].append(read)
        elif seq_method == "single_end":
            out["single"].append(read)

    return out

def get_seq_method(subsample):
    seq_method = get_subsample_attributes(subsample, "seq_method", pep)
    if isinstance(seq_method, list):
        seq_method = seq_method[0]
    return seq_method


# Symlink reference genome
rule symlink_genome:
    input: config["genome"]["path"]
    output: f"resources/symlink_genome/{config['genome']['name']}.fa"
    log: "logs/resources/symlink_genome.log"
    shell:
        """
        ln -s {input} {output} 2> {log}
        """

# Symlink reference genome annotations in gff3 format
rule symlink_gff3:
    input: config["gff3"]["path"]
    output: f"resources/symlink_gff3/{config['genome']['name']}.gff3"
    log: "logs/resources/symlink_gff3.log"
    shell:
        """
        ln -s {input} {output} 2> {log}
        """

# Create symlink that has the file renamed to specific ending
rule symlink_files:
    input:
        read1=lambda wildcards: get_fastqs(wildcards).get("r1", None),
        read2=lambda wildcards: get_fastqs(wildcards).get("r2", None)
    output:
        read1="resources/{project}/raw/{subsample}_R1.fastq.gz",
        read2="resources/{project}/raw/{subsample}_R2.fastq.gz"
    params:
        outdir="resources/{project}/raw/"
    log: "logs/{project}/raw/{subsample}.log"
    resources:
        mem_mb=1024
    run:
        if input.read1 and input.read2:
            # Create symlinks for paired-end reads if they exist
            shell("ln -s {input.read1} {output.read1} 2> {log}")
            shell("ln -s {input.read2} {output.read2} 2>> {log}")
        else:
            # If paired-end reads are missing, raise an error to skip the rule
            raise ValueError(f"Skipping symlink_files: Paired-end reads not found for {wildcards.subsample}.")

rule symlink_files_single:
    input: lambda wildcards: get_fastqs(wildcards)["single"],
    output: "resources/{project}/raw/{subsample}.fastq.gz",
    params:
        outdir="resources/{project}/raw/"
    log: "logs/{project}/raw/{subsample}.log"
    shell:
        """
        ln -s {input} {output} 2> {log}
        """
