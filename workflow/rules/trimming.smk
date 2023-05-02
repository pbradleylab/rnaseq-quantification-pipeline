""" Add rules to this section that are related to trimming.
"""

def getPairedFastqs(wildcards):
    out = {} # Dictionary that will hold two reads with r1 in index 0 and r2 in index 1
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    r1=[x for x in reads if ("_R1" in x or ".R1" in x or ".r1" in x or "_r1" in x or "_2.fq" in x or "_1.fastq" in x)]
    r2=[x for x in reads if ("_R2" in x or ".R2" in x or ".r2" in x or "_r2" in x or "_2.fq" in x or "_2.fastq" in x)]
    out["r1"] = r1[0]; out["r2"] = r2[0]
    return out

# rim Galore is a wrapper around Cutadapt and FastQC 
# to consistently apply adapter and quality trimming to 
# FastQ files, with extra functionality for RRBS data.
# Official documentation can be found here: https:/github.com/FelixKrueger/TrimGalore
rule trim_galore:
    input: unpack(getPairedFastqs)
    output:
        r1="results/trimming/trim_galore/{subsample}_R1_val_1.fq.gz",
        r2="results/trimming/trim_galore/{subsample}_R2_val_2.fq.gz"
    params:
        out_dir="results/trimming/trim_galore/" 
    log: "logs/trimming/trim_galore/{subsample}.log"
    resources:mem_mb=config["trim_galore"]["mem"]
    conda: "../envs/trimming.yml"
    shell: "trim_galore --paired {input} -o {params.out_dir}"