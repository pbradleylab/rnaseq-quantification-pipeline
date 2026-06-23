#!/usr/bin/env python3
"""Build a per-sample QC summary table from workflow outputs."""

import argparse
import csv
import re
import zipfile
from pathlib import Path


STATUS_RANK = {"PASS": 0, "WARN": 1, "FAIL": 2}


def sample_from_fastqc_filename(filename):
    sample = re.sub(r"\.(fastq|fq)(\.gz)?$", "", filename)
    sample = re.sub(r"(_R[12](_val_[12])?|_[12]_val_[12]|_trimmed)$", "", sample)
    return sample


def merge_status(current, new):
    if not current:
        return new
    return current if STATUS_RANK[current] >= STATUS_RANK[new] else new


def parse_fastqc(paths):
    raw = {}
    trimmed = {}
    for path in paths:
        path = Path(path)
        bucket = trimmed if "trim_galore" in str(path) or "_val_" in path.name or "_trimmed" in path.name else raw
        with zipfile.ZipFile(path) as archive:
            summary_names = [name for name in archive.namelist() if name.endswith("summary.txt")]
            for summary_name in summary_names:
                with archive.open(summary_name) as handle:
                    for raw_line in handle:
                        status, _module, filename = raw_line.decode().rstrip("\n").split("\t")
                        sample = sample_from_fastqc_filename(filename)
                        bucket[sample] = merge_status(bucket.get(sample), status)
    return raw, trimmed


def parse_fastq_screen(paths):
    out = {}
    for path in paths:
        path = Path(path)
        sample = path.parent.name
        out[sample] = "present"
    return out


def parse_star_logs(paths):
    out = {}
    for path in paths:
        path = Path(path)
        sample = path.parent.name
        metrics = {}
        with open(path) as handle:
            for line in handle:
                if "|" not in line:
                    continue
                key, value = [field.strip() for field in line.split("|", 1)]
                metrics[key] = value
        out[sample] = {
            "star_input_reads": metrics.get("Number of input reads", ""),
            "star_uniquely_mapped_percent": metrics.get("Uniquely mapped reads %", ""),
            "star_multimapped_percent": metrics.get("% of reads mapped to multiple loci", ""),
            "star_unmapped_too_short_percent": metrics.get("% of reads unmapped: too short", ""),
        }
    return out


def read_counts(path):
    counts = {}
    if not path:
        return counts
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        skip = {"target_id", "length"}
        for row in reader:
            for sample, value in row.items():
                if sample in skip or value in (None, ""):
                    continue
                counts[sample] = counts.get(sample, 0.0) + float(value)
    return counts


def read_metadata(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return [], []
    sample_col = "sample_name" if "sample_name" in rows[0] else "sample"
    return rows, [sample_col] + [col for col in rows[0] if col != sample_col]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fastqc", nargs="*", default=[])
    parser.add_argument("--fastq-screen", nargs="*", default=[])
    parser.add_argument("--star-logs", nargs="*", default=[])
    args = parser.parse_args()

    metadata_rows, metadata_cols = read_metadata(args.metadata)
    sample_col = metadata_cols[0]
    raw_fastqc, trimmed_fastqc = parse_fastqc(args.fastqc)
    fastq_screen = parse_fastq_screen(args.fastq_screen)
    star_logs = parse_star_logs(args.star_logs)
    assigned_counts = read_counts(args.counts)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = (
        metadata_cols
        + [
            "deseq2_included",
            "assigned_counts",
            "raw_fastqc_status",
            "trimmed_fastqc_status",
            "fastq_screen_report",
            "star_input_reads",
            "star_uniquely_mapped_percent",
            "star_multimapped_percent",
            "star_unmapped_too_short_percent",
        ]
    )
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for metadata_row in metadata_rows:
            sample = metadata_row[sample_col]
            star = star_logs.get(sample, {})
            row = {col: metadata_row.get(col, "") for col in metadata_cols}
            row.update(
                {
                    "deseq2_included": str(sample in assigned_counts).lower(),
                    "assigned_counts": assigned_counts.get(sample, ""),
                    "raw_fastqc_status": raw_fastqc.get(sample, ""),
                    "trimmed_fastqc_status": trimmed_fastqc.get(sample, ""),
                    "fastq_screen_report": fastq_screen.get(sample, ""),
                    "star_input_reads": star.get("star_input_reads", ""),
                    "star_uniquely_mapped_percent": star.get("star_uniquely_mapped_percent", ""),
                    "star_multimapped_percent": star.get("star_multimapped_percent", ""),
                    "star_unmapped_too_short_percent": star.get("star_unmapped_too_short_percent", ""),
                }
            )
            writer.writerow(row)


if __name__ == "__main__":
    main()
