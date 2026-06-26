""" Add rules to this section that are general utility functions or specific utilities.
"""
from pathlib import Path


def create_symlink(src, dest, log):
    src = Path(str(src)).resolve(strict=True)
    dest = Path(str(dest))
    log = Path(str(log))

    dest.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("")

    if dest.is_symlink():
        dest.unlink()
    elif dest.exists():
        raise FileExistsError(f"Refusing to replace non-symlink output: {dest}")

    dest.symlink_to(src)


# Symlink reference genome
rule symlink_genome:
    input: config["genome"]["path"]
    output: f"resources/symlink_genome/{config['genome']['name']}.fa"
    log: "logs/resources/symlink_genome.log"
    run:
        create_symlink(input[0], output[0], log[0])

# Symlink reference genome annotations in gff3 format
rule symlink_gff3:
    input: config["gff3"]["path"]
    output: f"resources/symlink_gff3/{config['genome']['name']}.gff3"
    log: "logs/resources/symlink_gff3.log"
    run:
        create_symlink(input[0], output[0], log[0])

# Create symlink that has the file renamed to specific ending
rule symlink_files:
    input:
        read1=lambda wildcards: get_fastqs(wildcards.subsample, pep).get("r1", None),
        read2=lambda wildcards: get_fastqs(wildcards.subsample, pep).get("r2", None)
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
            create_symlink(input.read1, output.read1, log[0])
            create_symlink(input.read2, output.read2, log[0])
        else:
            # If paired-end reads are missing, raise an error to skip the rule
            raise ValueError(f"Skipping symlink_files: Paired-end reads not found for {wildcards.subsample}.")

rule symlink_files_single:
    input: lambda wildcards: get_fastqs(wildcards.subsample, pep)["single"],
    output: "resources/{project}/raw/{subsample}.fastq.gz",
    wildcard_constraints:
        subsample=".+(?<!_R1)(?<!_R2)",
    params:
        outdir="resources/{project}/raw/"
    log: "logs/{project}/raw/{subsample}.log"
    run:
        create_symlink(input[0], output[0], log[0])
