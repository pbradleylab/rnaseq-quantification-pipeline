""" Add rules to this section that are related to trimming.
"""
# Trim Galore is a wrapper around Cutadapt and FastQC 
# to consistently apply adapter and quality trimming to 
# FastQ files, with extra functionality for RRBS data.
# Official documentation can be found here: https:/github.com/FelixKrueger/TrimGalore
rule trim_galore:
    input: 
        read1=rules.symlink_files.output[0],
        read2=rules.symlink_files.output[1]
    output:
        r1="results/{project}/trimming/trim_galore/{subsample}_R1_val_1.fq.gz",
        r2="results/{project}/trimming/trim_galore/{subsample}_R2_val_2.fq.gz"
    params:
        outdir="results/{project}/trimming/trim_galore/" 
    log: "logs/{project}/trimming/trim_galore/{subsample}.log"
    resources:mem_mb=config["trim_galore"]["mem"]
    conda: "../envs/trimming.yml"
    shell: "trim_galore --paired {input.read1} {input.read2} -o {params.outdir} 2> {log}"