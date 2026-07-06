import csv
import re
import sys
from pathlib import Path


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def as_list(value):
    return _as_list(value)


def config_flag(config, key, default="false"):
    return str(config.get(key, default)).lower() == "true"


def _format_list(values):
    return ", ".join(str(value) for value in sorted(values))


def _duplicates(values):
    seen = set()
    duplicates = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return duplicates


def read_deseq_metadata(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = reader.fieldnames or []
    return rows, columns


def get_design_formula_terms(design_formula):
    formula = str(design_formula).strip()
    if not formula.startswith("~"):
        return formula, []
    expression = formula[1:]
    expression = re.sub(r"`([^`]+)`", r"\1", expression)
    tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", expression)
    reserved = {"I", "factor", "as.factor", "log", "log2", "log10", "sqrt"}
    return formula, [token for token in tokens if token not in reserved]


def warn_preflight(message):
    print(f"WARNING: {message}", file=sys.stderr)


def validate_preflight_inputs(pep, config):
    errors = []

    metadata_path = config.get("metadata")
    if not metadata_path:
        errors.append("config/config.json must define a metadata CSV path.")
        metadata_rows = []
        metadata_columns = []
    elif not Path(metadata_path).exists():
        errors.append(f"Metadata file does not exist: {metadata_path}")
        metadata_rows = []
        metadata_columns = []
    else:
        metadata_rows, metadata_columns = read_deseq_metadata(metadata_path)
        if not metadata_rows:
            errors.append(f"Metadata file has no sample rows: {metadata_path}")

    sample_id_col = None
    if "sample_name" in metadata_columns:
        sample_id_col = "sample_name"
    elif "sample" in metadata_columns:
        sample_id_col = "sample"
    elif metadata_columns:
        errors.append("Metadata must include a sample_name or sample column.")

    metadata_samples = []
    if sample_id_col:
        metadata_samples = [
            row.get(sample_id_col, "").strip()
            for row in metadata_rows
            if row.get(sample_id_col, "").strip()
        ]
        duplicate_metadata_samples = _duplicates(metadata_samples)
        if duplicate_metadata_samples:
            errors.append(
                "Metadata contains duplicate sample IDs: "
                f"{_format_list(duplicate_metadata_samples)}"
            )

    deseq_config = config.get("deseq2", {})
    variable_to_analyze = deseq_config.get("variable_to_analyze")
    reference_in_variable = deseq_config.get("reference_in_variable")
    design_formula = deseq_config.get("design_formula")
    transform_method = str(deseq_config.get("transform_method", "vst")).lower()
    if transform_method not in {"vst", "rlog", "auto"}:
        errors.append("deseq2.transform_method must be one of: vst, rlog, auto.")

    strandedness_check = config.get("strandedness_check", {})
    strandedness_check_mode = str(strandedness_check.get("mode", "warn")).lower()
    if strandedness_check_mode not in {"warn", "fail"}:
        errors.append("strandedness_check.mode must be one of: warn, fail.")
    try:
        strandedness_min_fraction = float(
            strandedness_check.get("min_fraction", 0.75)
        )
    except (TypeError, ValueError):
        strandedness_min_fraction = 0.75
        errors.append("strandedness_check.min_fraction must be numeric.")
    if not 0 <= strandedness_min_fraction <= 1:
        errors.append("strandedness_check.min_fraction must be between 0 and 1.")
    try:
        strandedness_min_ratio = float(strandedness_check.get("min_ratio", 2.0))
    except (TypeError, ValueError):
        strandedness_min_ratio = 2.0
        errors.append("strandedness_check.min_ratio must be numeric.")
    if strandedness_min_ratio < 1:
        errors.append("strandedness_check.min_ratio must be at least 1.")

    sample_identity = config.get("sample_identity", {})
    sample_identity_group_column = sample_identity.get("group_column", "")
    if (
        sample_identity_group_column
        and metadata_columns
        and sample_identity_group_column not in metadata_columns
    ):
        errors.append(
            "sample_identity.group_column is not present in metadata: "
            f"{sample_identity_group_column}"
        )
    try:
        sample_identity_top_features = int(
            sample_identity.get("top_variable_features", 5000)
        )
    except (TypeError, ValueError):
        sample_identity_top_features = 5000
        errors.append("sample_identity.top_variable_features must be an integer.")
    if sample_identity_top_features < 1:
        errors.append("sample_identity.top_variable_features must be at least 1.")
    try:
        sample_identity_min_corr = float(
            sample_identity.get("min_nearest_correlation", 0.90)
        )
    except (TypeError, ValueError):
        sample_identity_min_corr = 0.90
        errors.append("sample_identity.min_nearest_correlation must be numeric.")
    if not -1 <= sample_identity_min_corr <= 1:
        errors.append(
            "sample_identity.min_nearest_correlation must be between -1 and 1."
        )
    try:
        sample_identity_same_group_margin = float(
            sample_identity.get("same_group_margin", 0.03)
        )
    except (TypeError, ValueError):
        sample_identity_same_group_margin = 0.03
        errors.append("sample_identity.same_group_margin must be numeric.")
    if sample_identity_same_group_margin < 0:
        errors.append("sample_identity.same_group_margin must be at least 0.")
    try:
        sample_identity_duplicate_min_corr = float(
            sample_identity.get("duplicate_min_correlation", 0.995)
        )
    except (TypeError, ValueError):
        sample_identity_duplicate_min_corr = 0.995
        errors.append("sample_identity.duplicate_min_correlation must be numeric.")
    if not -1 <= sample_identity_duplicate_min_corr <= 1:
        errors.append(
            "sample_identity.duplicate_min_correlation must be between -1 and 1."
        )
    try:
        sample_identity_duplicate_max_library_diff = float(
            sample_identity.get("duplicate_max_library_size_difference", 0.05)
        )
    except (TypeError, ValueError):
        sample_identity_duplicate_max_library_diff = 0.05
        errors.append(
            "sample_identity.duplicate_max_library_size_difference must be numeric."
        )
    if not 0 <= sample_identity_duplicate_max_library_diff <= 1:
        errors.append(
            "sample_identity.duplicate_max_library_size_difference must be "
            "between 0 and 1."
        )

    star_strandedness = str(
        config.get("star", {}).get("gene_counts_strandedness", "unstranded")
    ).lower()
    if star_strandedness not in {"unstranded", "forward", "reverse"}:
        errors.append(
            "star.gene_counts_strandedness must be one of: unstranded, "
            "forward, reverse."
        )
    featurecounts_strandedness = config.get("featurecounts", {}).get(
        "strandedness", 0
    )
    if str(featurecounts_strandedness) not in {"0", "1", "2"}:
        errors.append("featurecounts.strandedness must be one of: 0, 1, 2.")

    if config_flag(config, "run_sortmerna"):
        sortmerna_refs = _as_list(config.get("sortmerna", {}).get("refs", []))
        if not sortmerna_refs:
            errors.append(
                "run_sortmerna is true, but sortmerna.refs does not list any "
                "local rRNA reference FASTA files."
            )
        missing_sortmerna_refs = [
            ref for ref in sortmerna_refs if not Path(ref).exists()
        ]
        if missing_sortmerna_refs:
            errors.append(
                "SortMeRNA rRNA reference FASTA paths do not exist: "
                f"{_format_list(missing_sortmerna_refs)}"
            )

    if not design_formula and variable_to_analyze:
        design_formula = f"~ {variable_to_analyze}"
    if not design_formula:
        errors.append("config/tools.json must define deseq2.design_formula.")
    else:
        design_formula, design_terms = get_design_formula_terms(design_formula)
        if not design_formula.startswith("~"):
            errors.append("deseq2.design_formula must start with '~'.")
        missing_design_terms = [
            term for term in design_terms if term not in metadata_columns
        ]
        if missing_design_terms:
            errors.append(
                "DESeq2 design formula contains metadata columns that do not exist: "
                f"{_format_list(missing_design_terms)}"
            )
        if variable_to_analyze and variable_to_analyze not in design_terms:
            errors.append(
                f"DESeq2 design formula '{design_formula}' must include "
                f"variable_to_analyze '{variable_to_analyze}'."
            )
    try:
        min_replicates = int(deseq_config.get("min_replicates_per_condition", 2))
    except (TypeError, ValueError):
        min_replicates = 2
        errors.append("deseq2.min_replicates_per_condition must be an integer.")
    if min_replicates < 1:
        errors.append("deseq2.min_replicates_per_condition must be at least 1.")
    if not variable_to_analyze:
        errors.append("config/tools.json must define deseq2.variable_to_analyze.")
    elif metadata_columns and variable_to_analyze not in metadata_columns:
        errors.append(
            f"DESeq2 condition column '{variable_to_analyze}' is missing from metadata."
        )

    if variable_to_analyze in metadata_columns:
        level_counts = {}
        for row in metadata_rows:
            level = row.get(variable_to_analyze, "").strip()
            sample = row.get(sample_id_col, "").strip() if sample_id_col else ""
            if level and sample:
                level_counts[level] = level_counts.get(level, 0) + 1
        levels = set(level_counts)
        if len(levels) < 2:
            errors.append(
                f"DESeq2 condition column '{variable_to_analyze}' must contain "
                "at least two levels."
            )
        if not reference_in_variable:
            errors.append("config/tools.json must define deseq2.reference_in_variable.")
        elif reference_in_variable not in levels:
            errors.append(
                f"DESeq2 reference level '{reference_in_variable}' is not present in "
                f"metadata column '{variable_to_analyze}'. Observed levels: "
                f"{_format_list(levels)}"
            )
        low_replicate_levels = {
            level: count
            for level, count in level_counts.items()
            if count < min_replicates
        }
        if low_replicate_levels:
            level_summary = ", ".join(
                f"{level}={count}"
                for level, count in sorted(low_replicate_levels.items())
            )
            warn_preflight(
                f"DESeq2 condition column '{variable_to_analyze}' has fewer "
                f"than {min_replicates} biological replicates for at least one "
                f"condition ({level_summary}). DESeq2 may run, but interpretation "
                "will be weak."
            )

    pep_sample_names = pep.sample_table.sample_name.tolist()
    duplicate_pep_samples = _duplicates(pep_sample_names)
    if duplicate_pep_samples:
        errors.append(
            "PEP sample table contains duplicate sample IDs: "
            f"{_format_list(duplicate_pep_samples)}"
        )

    pep_samples = set(pep_sample_names)
    subsample_samples = set(pep.subsample_table.sample_name.tolist())
    samples_without_subsamples = pep_samples - subsample_samples
    if samples_without_subsamples:
        errors.append(
            "PEP samples without subsample FASTQ entries: "
            f"{_format_list(samples_without_subsamples)}"
        )

    metadata_sample_set = set(metadata_samples)
    samples_without_metadata = subsample_samples - metadata_sample_set
    metadata_without_fastqs = metadata_sample_set - subsample_samples
    if samples_without_metadata:
        errors.append(
            "PEP FASTQ samples missing from metadata: "
            f"{_format_list(samples_without_metadata)}"
        )
    if metadata_without_fastqs:
        errors.append(
            "Metadata samples missing from PEP FASTQ samples: "
            f"{_format_list(metadata_without_fastqs)}"
        )

    for subsample in pep.subsample_table.subsample.tolist():
        subsample_rows = pep.subsample_table[
            pep.subsample_table["subsample"] == subsample
        ]
        sample_name = subsample_rows["sample_name"].tolist()[0]
        seq_method = get_seq_method(subsample, pep)
        fastqs = get_fastqs(subsample, pep)
        reads = _as_list(get_subsample_attributes(subsample, "reads", pep))

        missing_reads = [read for read in reads if not Path(read).exists()]
        if missing_reads:
            errors.append(
                f"Subsample '{subsample}' has FASTQ paths that do not exist: "
                f"{_format_list(missing_reads)}"
            )

        if seq_method == "paired_end":
            if not fastqs["r1"] or not fastqs["r2"]:
                errors.append(
                    f"Paired-end subsample '{subsample}' for sample '{sample_name}' "
                    "must have both R1 and R2 FASTQs."
                )
            elif len(fastqs["r1"]) != len(fastqs["r2"]):
                errors.append(
                    f"Paired-end subsample '{subsample}' for sample '{sample_name}' "
                    "has a different number of R1 and R2 FASTQs."
                )
        elif seq_method == "single_end":
            if not fastqs["single"]:
                errors.append(
                    f"Single-end subsample '{subsample}' for sample '{sample_name}' "
                    "must have at least one FASTQ."
                )
        else:
            errors.append(
                f"Subsample '{subsample}' has unsupported seq_method '{seq_method}'."
            )

    if errors:
        raise ValueError(
            "Preflight input validation failed:\n- " + "\n- ".join(errors)
        )


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
            out.append(rules.kallisto.output.h5.format(project=project, subsample=subsample))
        elif seq_method == "single_end":
            out.append(rules.kallisto_single.output.h5.format(project=project, subsample=subsample))
    return out


def get_quantification_reads(wildcards, rules, config):
    if config_flag(config, "run_sortmerna"):
        return [
            rules.sortmerna.output.r1.format(
                project=wildcards.project, subsample=wildcards.subsample
            ),
            rules.sortmerna.output.r2.format(
                project=wildcards.project, subsample=wildcards.subsample
            ),
        ]
    return [
        rules.trim_galore.output.r1.format(
            project=wildcards.project, subsample=wildcards.subsample
        ),
        rules.trim_galore.output.r2.format(
            project=wildcards.project, subsample=wildcards.subsample
        ),
    ]


def get_quantification_reads_single(wildcards, rules, config):
    if config_flag(config, "run_sortmerna"):
        return rules.sortmerna_single.output.reads.format(
            project=wildcards.project, subsample=wildcards.subsample
        )
    return rules.trim_galore_single.output[0].format(
        project=wildcards.project, subsample=wildcards.subsample
    )


def get_kallisto_abundance_tsv(wildcards, pep, rules):
    out = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        if project != wildcards.project:
            continue
        seq_method = get_seq_method(subsample, pep)
        if seq_method == "paired_end":
            out.append(rules.kallisto.output.tsv.format(project=project, subsample=subsample))
        elif seq_method == "single_end":
            out.append(rules.kallisto_single.output.tsv.format(project=project, subsample=subsample))
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


def get_strandedness_reports(wildcards, pep, rules, config):
    out = []
    if not config_flag(config.get("strandedness_check", {}), "enabled", "true"):
        return out
    if config["quantification_tool"].lower() not in ["star", "featurecounts"]:
        return out
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        if project != wildcards.project:
            continue
        seq_method = get_seq_method(subsample, pep)
        if seq_method == "paired_end":
            out.append(
                rules.strandedness_check.output[0].format(
                    project=project, subsample=subsample
                )
            )
        elif seq_method == "single_end":
            out.append(
                rules.strandedness_check_single.output[0].format(
                    project=project, subsample=subsample
                )
            )
    return out


def get_featurecounts_counts(wildcards, pep, rules):
    out = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        if project != wildcards.project:
            continue
        seq_method = get_seq_method(subsample, pep)
        if seq_method == "paired_end":
            out.append(rules.featurecounts.output.counts.format(project=project, subsample=subsample))
        elif seq_method == "single_end":
            out.append(rules.featurecounts_single.output.counts.format(project=project, subsample=subsample))
    return out
