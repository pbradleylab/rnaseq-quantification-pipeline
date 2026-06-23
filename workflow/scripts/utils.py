def get_subsample_attributes(subsample, attribute, pep):
    subsample_rows = pep.subsample_table[pep.subsample_table["subsample"] == subsample]
    return pep.get_sample(subsample_rows["sample_name"].tolist()[0])[attribute]


def get_seq_method(subsample, pep):
    seq_method = get_subsample_attributes(subsample, "seq_method", pep)
    if isinstance(seq_method, list):
        seq_method = seq_method[0]
    return seq_method


def get_fastqs(subsample, pep):
    out = {"r1": [], "r2": [], "single": []}
    reads = get_subsample_attributes(subsample, "reads", pep)
    if isinstance(reads, str):
        reads = [reads]

    seq_method = get_seq_method(subsample, pep)
    for read in reads:
        if seq_method == "paired_end":
            if any(x in read for x in ["_R1", ".R1", ".r1", "_r1", "_1.fq", "_1.fastq"]):
                out["r1"].append(read)
            elif any(x in read for x in ["_R2", ".R2", ".r2", "_r2", "_2.fq", "_2.fastq"]):
                out["r2"].append(read)
        elif seq_method == "single_end":
            out["single"].append(read)

    return out


def get_kallisto_h5(wildcards, pep, rules):
    out = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        if project != wildcards.project:
            continue
        seq_method = get_seq_method(subsample, pep)
        if seq_method == "paired_end":
            out.append(rules.kallisto.output[0].format(project=project, subsample=subsample))
        elif seq_method == "single_end":
            out.append(rules.kallisto_single.output[0].format(project=project, subsample=subsample))
    return out


def get_star_gene_counts(wildcards, pep, rules):
    out = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        if project != wildcards.project:
            continue
        seq_method = get_seq_method(subsample, pep)
        if seq_method == "paired_end":
            out.append(rules.star_reads_per_gene.output[0].format(project=project, subsample=subsample))
        elif seq_method == "single_end":
            out.append(rules.star_reads_per_gene_single.output[0].format(project=project, subsample=subsample))
    return out
