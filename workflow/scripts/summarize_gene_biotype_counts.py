#!/usr/bin/env python3
"""Summarize count support by annotation gene_biotype and count class."""

import argparse
import csv
import sys
from collections import defaultdict


COUNT_CLASSES = (
    ("zero", 0, 0),
    ("low_1_10", 1, 10),
    ("moderate_11_100", 11, 100),
    ("high_gt_100", 101, None),
)


def parse_gff_attributes(attribute_text):
    attributes = {}
    for field in attribute_text.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attributes[key.strip()] = value.strip()
    return attributes


def read_annotation(path, feature_type, id_attribute):
    annotations = {}
    saw_biotype = False
    with open(path, newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != feature_type:
                continue
            attributes = parse_gff_attributes(fields[8])
            feature_id = attributes.get(id_attribute)
            if not feature_id:
                continue
            biotype = attributes.get("gene_biotype")
            if biotype:
                saw_biotype = True
            annotations[feature_id] = {
                "gene_biotype": biotype or "",
            }
    return annotations, saw_biotype


def read_count_totals(path):
    totals = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "target_id" not in reader.fieldnames:
            raise SystemExit(f"Count matrix must contain a target_id column: {path}")
        sample_columns = [
            column for column in reader.fieldnames if column not in {"target_id", "length"}
        ]
        for row in reader:
            target_id = row.get("target_id", "")
            if not target_id:
                continue
            total = 0.0
            for sample in sample_columns:
                value = row.get(sample, "0") or "0"
                total += float(value)
            totals[target_id] = total
    return totals


def count_class(total):
    for label, lower, upper in COUNT_CLASSES:
        if upper is None and total >= lower:
            return label
        if upper is not None and lower <= total <= upper:
            return label
    return "unclassified"


def write_message(path, message):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["status", "message"])
        writer.writerow(["not_calculated", message])
    print(message, file=sys.stderr)


def write_summary(path, annotations, count_totals, feature_type):
    summary = defaultdict(lambda: {"genes": 0, "genes_in_count_matrix": 0, "total_counts": 0.0})
    for feature_id, annotation in annotations.items():
        total = count_totals.get(feature_id, 0.0)
        key = (annotation["gene_biotype"], count_class(total))
        summary[key]["genes"] += 1
        summary[key]["genes_in_count_matrix"] += int(feature_id in count_totals)
        summary[key]["total_counts"] += total

    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "gene_biotype",
                "count_class",
                "annotation_genes",
                "genes_in_count_matrix",
                "total_counts",
            ]
        )
        for (biotype, count_class_label), values in sorted(summary.items()):
            writer.writerow(
                [
                    biotype,
                    count_class_label,
                    values["genes"],
                    values["genes_in_count_matrix"],
                    f"{values['total_counts']:.6f}",
                ]
            )

    print(
        "Gene biotype/count class summary calculated for "
        f"{len(annotations)} annotated {feature_type} records.",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--feature-type", default="gene")
    parser.add_argument("--id-attribute", default="ID")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    annotations, saw_biotype = read_annotation(
        args.annotation, args.feature_type, args.id_attribute
    )
    if not annotations:
        raise SystemExit(
            "No annotation records found with "
            f"feature_type='{args.feature_type}' and id_attribute='{args.id_attribute}'."
        )
    if not saw_biotype:
        write_message(
            args.output,
            "Annotation does not contain gene_biotype attributes, so gene "
            "biotype/count class summaries were not calculated.",
        )
        return

    count_totals = read_count_totals(args.counts)
    write_summary(args.output, annotations, count_totals, args.feature_type)


if __name__ == "__main__":
    main()
