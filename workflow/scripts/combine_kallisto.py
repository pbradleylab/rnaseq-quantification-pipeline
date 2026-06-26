#!/usr/bin/env python3
"""Combine Kallisto abundance.tsv files into count and TPM matrices."""

import argparse
import csv
from pathlib import Path


def abundance_files(outdirs):
    files = []
    for outdir in outdirs:
        files.extend(sorted(Path(outdir).glob("*/abundance.tsv")))
    return files


def read_abundance(path):
    rows = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["target_id"]] = row
    return rows


def write_matrix(path, samples, abundances, value_col):
    targets = sorted(set().union(*(sample.keys() for sample in abundances)))
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target_id", "length", *samples])
        for target in targets:
            length = next((rows[target]["length"] for rows in abundances if target in rows), 0)
            values = [rows.get(target, {}).get(value_col, 0) for rows in abundances]
            writer.writerow([target, length, *values])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--tpm", required=True)
    parser.add_argument("outdirs", nargs="+")
    args = parser.parse_args()

    files = abundance_files(args.outdirs)
    if not files:
        raise SystemExit("No abundance.tsv files found.")

    samples = [path.parent.name for path in files]
    abundances = [read_abundance(path) for path in files]
    write_matrix(args.counts, samples, abundances, "est_counts")
    write_matrix(args.tpm, samples, abundances, "tpm")


if __name__ == "__main__":
    main()
