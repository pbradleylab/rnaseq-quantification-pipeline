from scripts.utils import *
configfile: "config/config.json"

def get_fastq_input(wildcards):
    return get_subsample_attributes(wildcards.subsample, "reads", pep)

def get_multiqc_subsamples(wildcards):
    metricsLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        
        # Always run rules on the outside
        metricsLST.append(rules.fastqc.output[0].format(project=project, subsample=subsample))
        metricsLST.append(rules.fastq_screen.output[0].format(project=project, subsample=subsample))
        metricsLST.append(rules.combine_quants.output[0].format(project=project))

    return metricsLST

def get_combine_quant_input(wildcards):
    metricsLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        metricsLST.append(rules.kallisto.output[0].format(project=project, subsample=subsample))

    return metricsLST

# Calculates the contamination accross in the samples.
# The config must be supplied when calling this rule.
rule fastq_screen:
    input: unpack(get_fastq_input)
    output:
        html="../metrics/{project}/fastqScreen/{subsample}_screen.html",
        txt="../metrics/{project}/fastqScreen/{subsample}_screen.txt",
    params:
        fastq_screen_config="../config/fastq_screen_config.conf",
        subset=config["fastq_screen"]["subset"],
        aligner=config["fastq_screen"]["aligner"],
        tempDir="../metrics/{project}/fastqScreen/",
        sample="{subsample}"
    log: "../logs/{project}/metrics/fastqScreen/{subsample}.log"
    threads: config["fastq_screen"]["threads"]
    resources: mem_mb = config["fastq_screen"]["mem"]
    conda: "../envs/quality_control.yml"
    shell: 
        """
        fastq_screen --outdir {params.tempDir} --force --aligner {params.aligner} --conf {params.fastq_screen_config} --subset {params.subset} --threads {threads} {input} &> {log}
        """

rule fastqc:
    input: unpack(get_fastq_input)
    output: "../metrics/{project}/fastqc/{subsample}/{subsample}.zip"
    log: "../logs/{project}/metrics/fastqc/{subsample}.log"
    params:
        out_dir="../metrics/{project}/fastqScreen/",
        subsample="{subsample}"
    threads: config["fastqc"]["threads"]
    resources: mem_mb = config["fastq_screen"]["mem"]
    conda: "../envs/quality_control.yml"
    shell: 
        """
        "fastqc {input} --extract -t {threads} -o {output} &> {log}"
        mv {output} {params.out_dir}/{params.subsample}.zip
        """

rule combine_quants:
    input: unpack(get_combine_quant_input)
    output: "../results/{project}/quantification/combineQuants/combine.quant.txt"
    params:
        ext=".kallisto_gene.txt",
        gtf=config["gff3"]
    resources: mem=config["combine_quants"]["mem"]
    run:
        dirs = list(set([os.path.dirname(x) for x in input]))
        file_dir = dirs[0]
        abuns = os.listdir(fileDir)
        abuns_pd = combine(abuns, file_dir, str(params.ext), str(params.gtf))
        # Write to a file
        out = open(str(output), 'w')
        abuns_pd.to_csv(out, index = False, sep='\t')

# Generates html reports of analysis / QC done using MultiQC
# as shown here: https://multiqc.info/
# Files that are too large may be skipped and this rule may
# fail even if the fastqc is successfully run. To adjust this
# the config file for fastqc may be adjusted.
rule multiqc:
    input:get_multiqc_subsamples
    output: 
        html="../metrics/{project}/multiqc/multiqc_report.html"    
    params:
        search=directory("../"),
        output_dir=directory("../metrics/{project}/multiqc/")
    log: "../logs/{project}/metrics/multiqc/multiqc_report.log"
    resources: mem_mb = config["multiqc"]["mem"]
    conda: "../envs/quality_control.yml"
    shell:
        """
        multiqc -f {params.search} --cl_config log_filesize_limit:2000000000 --outdir {params.output_dir} 2> {log}
        """