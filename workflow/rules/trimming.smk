
def getPairedFastqs(wildcards):
    out = {}
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    r1=[x for x in reads if ("_R1" in x or ".R1" in x or ".r1" in x or "_r1" in x or "_2.fq" in x or "_1.fastq" in x)]
    r2=[x for x in reads if ("_R2" in x or ".R2" in x or ".r2" in x or "_r2" in x or "_2.fq" in x or "_2.fastq" in x)]
    out["r1"] = r1[0]; out["r2"] = r2[0]
    return out

rule trim_galore:
    input: unpack(getPairedFastqs)
    output:
        r1="../results/trimming/trim_galore/{subsample}.R1.fastq.gz",
        r1_unpaired="../results/trimming/trim_galore/{subsample}.R1.unpaired.fastq.gz",
        r2="../results/trimming/trim_galore/{subsample}.R2.fastq.gz",
        r2_unpaired="../results/trimming/trim_galore/{subsample}.R2.unpaired.fastq.gz"        
    log: "../logs/trimming/trim_galore/{subsample}.log"
    resources:mem_mb=config["trim_galore"]["mem"]
    conda: "../envs/trimming.yml"
    shell: "trim_galore --paired *.clock_UMI.R1.fq.gz *.clock_UMI.R2.fq.gz"