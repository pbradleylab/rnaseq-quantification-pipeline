#!/usr/bin/env python3
"""Summarize sample identity from expression-profile similarity."""

import argparse
import csv
import hashlib
import math
from pathlib import Path


SKIP_COUNT_COLUMNS = {"length"}


def read_metadata(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    if not fieldnames:
        return [], [], ""
    sample_col = "sample_name" if "sample_name" in fieldnames else "sample"
    if sample_col not in fieldnames:
        sample_col = fieldnames[0]
    return rows, fieldnames, sample_col


def read_count_matrix(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return [], [], {}
        feature_col = reader.fieldnames[0]
        samples = [
            col
            for col in reader.fieldnames[1:]
            if col and col not in SKIP_COUNT_COLUMNS
        ]
        values = {sample: [] for sample in samples}
        features = []
        for row in reader:
            feature = row.get(feature_col, "")
            if not feature:
                continue
            features.append(feature)
            for sample in samples:
                try:
                    value = float(row.get(sample, 0) or 0)
                except ValueError:
                    value = 0.0
                values[sample].append(value)
    return features, samples, values


def variance(values):
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return sum((value - mean) ** 2 for value in values) / (len(values) - 1)


def pearson(left, right):
    if len(left) != len(right) or len(left) < 2:
        return None
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = 0.0
    left_ss = 0.0
    right_ss = 0.0
    for left_value, right_value in zip(left, right):
        left_delta = left_value - left_mean
        right_delta = right_value - right_mean
        numerator += left_delta * right_delta
        left_ss += left_delta**2
        right_ss += right_delta**2
    denominator = math.sqrt(left_ss * right_ss)
    if denominator == 0:
        return None
    return numerator / denominator


def select_variable_feature_indexes(features, samples, values, top_n):
    if not features or not samples:
        return []
    log_values = {
        sample: [math.log2(max(value, 0.0) + 1.0) for value in values[sample]]
        for sample in samples
    }
    ranked = []
    for idx, feature in enumerate(features):
        feature_values = [log_values[sample][idx] for sample in samples]
        ranked.append((variance(feature_values), feature, idx))
    ranked.sort(key=lambda item: (-item[0], item[1]))
    selected = [idx for feature_variance, _feature, idx in ranked if feature_variance > 0]
    if not selected:
        selected = [idx for _feature_variance, _feature, idx in ranked]
    if top_n > 0:
        selected = selected[:top_n]
    return selected


def selected_log_profiles(samples, values, indexes):
    profiles = {}
    for sample in samples:
        profiles[sample] = [
            math.log2(max(values[sample][idx], 0.0) + 1.0) for idx in indexes
        ]
    return profiles


def fingerprint(features, values, indexes, sample, top_n=25):
    ranked = sorted(
        (
            (values[sample][idx], features[idx])
            for idx in indexes
            if idx < len(values[sample])
        ),
        key=lambda item: (-item[0], item[1]),
    )[:top_n]
    payload = "|".join(f"{feature}:{value:.6g}" for value, feature in ranked)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def library_size_relative_difference(left, right):
    denominator = max(left, right)
    if denominator <= 0:
        return None
    return abs(left - right) / denominator


def metadata_differences(sample, other, metadata_by_sample, metadata_cols, sample_col):
    left = metadata_by_sample.get(sample, {})
    right = metadata_by_sample.get(other, {})
    shared = []
    different = []
    for col in metadata_cols:
        if col == sample_col:
            continue
        left_value = left.get(col, "")
        right_value = right.get(col, "")
        if left_value == "" or right_value == "":
            continue
        if left_value == right_value:
            shared.append(col)
        else:
            different.append(col)
    return shared, different


def nearest_neighbors(samples, profiles):
    correlations = {sample: {} for sample in samples}
    for i, sample in enumerate(samples):
        for other in samples[i + 1 :]:
            corr = pearson(profiles[sample], profiles[other])
            correlations[sample][other] = corr
            correlations[other][sample] = corr
    nearest = {}
    for sample in samples:
        valid = [
            (other, corr)
            for other, corr in correlations[sample].items()
            if corr is not None
        ]
        valid.sort(key=lambda item: (-item[1], item[0]))
        nearest[sample] = valid[0] if valid else ("", None)
    return nearest, correlations


def best_group_neighbor(sample, correlations, metadata_by_sample, group_col, same_group):
    if not group_col:
        return "", None
    group = metadata_by_sample.get(sample, {}).get(group_col, "")
    if not group:
        return "", None
    candidates = []
    for other, corr in correlations.get(sample, {}).items():
        if corr is None:
            continue
        other_group = metadata_by_sample.get(other, {}).get(group_col, "")
        if not other_group:
            continue
        if (other_group == group) == same_group:
            candidates.append((other, corr))
    candidates.sort(key=lambda item: (-item[1], item[0]))
    return candidates[0] if candidates else ("", None)


def fmt_corr(value):
    return "" if value is None else f"{value:.4f}"


def write_similarity_matrix(path, ordered_samples, count_samples, correlations):
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    count_sample_set = set(count_samples)
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", *ordered_samples])
        for sample in ordered_samples:
            row = [sample]
            for other in ordered_samples:
                if sample == other and sample in count_sample_set:
                    row.append("1.0000")
                elif sample in count_sample_set and other in count_sample_set:
                    row.append(fmt_corr(correlations.get(sample, {}).get(other)))
                else:
                    row.append("")
            writer.writerow(row)


def duplicate_library_candidates(
    sample,
    count_samples,
    correlations,
    library_totals,
    fingerprints,
    min_correlation,
    max_library_size_difference,
):
    if sample not in count_samples:
        return []
    candidates = []
    for other in count_samples:
        if other == sample:
            continue
        corr = correlations.get(sample, {}).get(other)
        if corr is None or corr < min_correlation:
            continue
        library_diff = library_size_relative_difference(
            library_totals.get(sample, 0.0), library_totals.get(other, 0.0)
        )
        if library_diff is None or library_diff > max_library_size_difference:
            continue
        fingerprint_match = fingerprints.get(sample) == fingerprints.get(other)
        candidates.append(
            {
                "sample": other,
                "correlation": corr,
                "library_size_relative_difference": library_diff,
                "fingerprint_match": fingerprint_match,
            }
        )
    candidates.sort(
        key=lambda item: (
            -item["correlation"],
            item["library_size_relative_difference"],
            item["sample"],
        )
    )
    return candidates


def identity_status(
    sample,
    count_samples,
    library_total,
    nearest_corr,
    same_group_corr,
    different_group_corr,
    min_nearest_correlation,
    same_group_margin,
):
    if sample not in count_samples:
        return "fail", "Sample is present in metadata but missing from count matrix."
    if library_total <= 0:
        return "fail", "Sample has zero assigned counts."
    if nearest_corr is None:
        return "limited", "No comparable expression neighbor could be calculated."
    if nearest_corr < min_nearest_correlation:
        return "review", "Nearest expression neighbor has low correlation."
    if (
        same_group_corr is not None
        and different_group_corr is not None
        and different_group_corr - same_group_corr >= same_group_margin
    ):
        return (
            "review",
            "Sample is closer to a different metadata group than to its "
            "closest same-group sample.",
        )
    return "pass", ""


def main():
    parser = argparse.ArgumentParser(
        description="Build an expression-profile sample identity report."
    )
    parser.add_argument("--counts", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--similarity-matrix", required=True)
    parser.add_argument("--group-column", default="")
    parser.add_argument("--top-variable-features", type=int, default=5000)
    parser.add_argument("--min-nearest-correlation", type=float, default=0.90)
    parser.add_argument("--same-group-margin", type=float, default=0.03)
    parser.add_argument("--duplicate-min-correlation", type=float, default=0.995)
    parser.add_argument(
        "--duplicate-max-library-size-difference", type=float, default=0.05
    )
    args = parser.parse_args()

    metadata_rows, metadata_cols, sample_col = read_metadata(args.metadata)
    features, count_samples, count_values = read_count_matrix(args.counts)
    metadata_by_sample = {
        row.get(sample_col, ""): row
        for row in metadata_rows
        if row.get(sample_col, "")
    }

    group_col = args.group_column if args.group_column in metadata_cols else ""
    ordered_samples = [row.get(sample_col, "") for row in metadata_rows if row.get(sample_col, "")]
    ordered_samples.extend(sample for sample in count_samples if sample not in metadata_by_sample)

    indexes = select_variable_feature_indexes(
        features, count_samples, count_values, args.top_variable_features
    )
    profiles = selected_log_profiles(count_samples, count_values, indexes)
    nearest, correlations = nearest_neighbors(count_samples, profiles)
    library_totals = {
        sample: sum(count_values.get(sample, [])) for sample in count_samples
    }
    fingerprints = {
        sample: fingerprint(features, count_values, indexes, sample)
        for sample in count_samples
    }
    write_similarity_matrix(
        args.similarity_matrix, ordered_samples, count_samples, correlations
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample",
        "metadata_present",
        "counts_present",
        "assigned_counts",
        "fingerprint_sha256_16",
        "features_used",
        "nearest_sample",
        "nearest_correlation",
        "nearest_shared_metadata_columns",
        "nearest_different_metadata_columns",
        "group_column",
        "group_value",
        "nearest_same_group_sample",
        "nearest_same_group_correlation",
        "nearest_different_group_sample",
        "nearest_different_group_correlation",
        "duplicate_library_candidate",
        "duplicate_candidate_count",
        "duplicate_candidate_samples",
        "duplicate_max_correlation",
        "duplicate_min_library_size_relative_difference",
        "duplicate_exact_fingerprint_match",
        "identity_status",
        "message",
    ]
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sample in ordered_samples:
            nearest_sample, nearest_corr = nearest.get(sample, ("", None))
            same_sample, same_corr = best_group_neighbor(
                sample, correlations, metadata_by_sample, group_col, same_group=True
            )
            different_sample, different_corr = best_group_neighbor(
                sample, correlations, metadata_by_sample, group_col, same_group=False
            )
            library_total = library_totals.get(sample, 0.0)
            status, message = identity_status(
                sample,
                set(count_samples),
                library_total,
                nearest_corr,
                same_corr,
                different_corr,
                args.min_nearest_correlation,
                args.same_group_margin,
            )
            if sample in count_samples and sample not in metadata_by_sample:
                status = "review"
                message = "Sample is present in count matrix but missing from metadata."
            duplicate_candidates = duplicate_library_candidates(
                sample,
                count_samples,
                correlations,
                library_totals,
                fingerprints,
                args.duplicate_min_correlation,
                args.duplicate_max_library_size_difference,
            )
            duplicate_candidate = bool(duplicate_candidates)
            if duplicate_candidate:
                duplicate_message = (
                    "Possible duplicate library: expression profile and assigned "
                    "count total closely match another sample."
                )
                if status == "pass":
                    status = "review"
                    message = duplicate_message
                elif message:
                    message = f"{message} {duplicate_message}"
                else:
                    message = duplicate_message
            shared, different = metadata_differences(
                sample, nearest_sample, metadata_by_sample, metadata_cols, sample_col
            )
            duplicate_sample_details = []
            for candidate in duplicate_candidates:
                library_diff = candidate["library_size_relative_difference"]
                duplicate_sample_details.append(
                    (
                        "{sample}:r={corr}:library_diff={library_diff}:"
                        "fingerprint_match={match}"
                    ).format(
                        sample=candidate["sample"],
                        corr=fmt_corr(candidate["correlation"]),
                        library_diff=f"{library_diff:.4f}",
                        match=str(candidate["fingerprint_match"]).lower(),
                    )
                )
            duplicate_min_library_diff = (
                min(
                    candidate["library_size_relative_difference"]
                    for candidate in duplicate_candidates
                )
                if duplicate_candidates
                else None
            )
            duplicate_exact_fingerprint_match = any(
                candidate["fingerprint_match"] for candidate in duplicate_candidates
            )
            writer.writerow(
                {
                    "sample": sample,
                    "metadata_present": str(sample in metadata_by_sample).lower(),
                    "counts_present": str(sample in count_samples).lower(),
                    "assigned_counts": (
                        f"{library_total:.6g}" if sample in count_samples else ""
                    ),
                    "fingerprint_sha256_16": (
                        fingerprints[sample]
                        if sample in count_samples
                        else ""
                    ),
                    "features_used": len(indexes),
                    "nearest_sample": nearest_sample,
                    "nearest_correlation": fmt_corr(nearest_corr),
                    "nearest_shared_metadata_columns": ",".join(shared),
                    "nearest_different_metadata_columns": ",".join(different),
                    "group_column": group_col,
                    "group_value": metadata_by_sample.get(sample, {}).get(group_col, ""),
                    "nearest_same_group_sample": same_sample,
                    "nearest_same_group_correlation": fmt_corr(same_corr),
                    "nearest_different_group_sample": different_sample,
                    "nearest_different_group_correlation": fmt_corr(different_corr),
                    "duplicate_library_candidate": str(duplicate_candidate).lower(),
                    "duplicate_candidate_count": len(duplicate_candidates),
                    "duplicate_candidate_samples": ";".join(duplicate_sample_details),
                    "duplicate_max_correlation": (
                        fmt_corr(duplicate_candidates[0]["correlation"])
                        if duplicate_candidates
                        else ""
                    ),
                    "duplicate_min_library_size_relative_difference": (
                        f"{duplicate_min_library_diff:.4f}"
                        if duplicate_min_library_diff is not None
                        else ""
                    ),
                    "duplicate_exact_fingerprint_match": str(
                        duplicate_exact_fingerprint_match
                    ).lower(),
                    "identity_status": status,
                    "message": message,
                }
            )


if __name__ == "__main__":
    main()
