#!/usr/bin/env python3
"""Combine STAR ReadsPerGene.out.tab files into a DESeq2 count matrix."""

import argparse
import csv
from pathlib import Path

COUNT_COLUMNS = {
    "unstranded": 1,
    "forward": 2,
    "reverse": 3,
}


def sample_name(path):
    return Path(path).parent.name


def read_counts(path, column):
    counts = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("N_"):
                continue
            counts[row[0]] = row[column]
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--strandedness",
        choices=sorted(COUNT_COLUMNS),
        default="unstranded",
        help=(
            "STAR ReadsPerGene.out.tab count column to use: unstranded, "
            "forward, or reverse."
        ),
    )
    parser.add_argument("reads_per_gene", nargs="+")
    args = parser.parse_args()

    column = COUNT_COLUMNS[args.strandedness]
    samples = [sample_name(path) for path in args.reads_per_gene]
    sample_counts = [read_counts(path, column) for path in args.reads_per_gene]
    genes = sorted(set().union(*(counts.keys() for counts in sample_counts)))

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target_id", *samples])
        for gene in genes:
            writer.writerow([gene, *(counts.get(gene, 0) for counts in sample_counts)])


if __name__ == "__main__":
    main()
