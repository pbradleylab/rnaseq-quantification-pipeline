#!/usr/bin/env python3
"""Combine per-sample featureCounts outputs into a count matrix."""

import argparse
import csv
from pathlib import Path


def sample_name(path):
    return Path(path).parent.name


def read_counts(path):
    counts = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = None
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if header is None:
                header = row
                continue
            counts[row[0]] = row[-1]
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("featurecounts", nargs="+")
    args = parser.parse_args()

    samples = [sample_name(path) for path in args.featurecounts]
    sample_counts = [read_counts(path) for path in args.featurecounts]
    genes = sorted(set().union(*(counts.keys() for counts in sample_counts)))

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target_id", *samples])
        for gene in genes:
            writer.writerow([gene, *(counts.get(gene, 0) for counts in sample_counts)])


if __name__ == "__main__":
    main()
