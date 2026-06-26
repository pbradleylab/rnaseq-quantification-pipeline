#!/usr/bin/env python3
"""Report overlap between count-matrix IDs and annotation feature IDs."""

import argparse
import csv
import sys


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


def read_count_ids(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "target_id" not in reader.fieldnames:
            raise SystemExit(f"Count matrix must contain a target_id column: {path}")
        return {row["target_id"] for row in reader if row.get("target_id")}


def read_annotation_ids(path, feature_type, attribute):
    ids = set()
    with open(path, newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if feature_type and fields[2] != feature_type:
                continue
            attributes = parse_gff_attributes(fields[8])
            value = attributes.get(attribute)
            if value:
                ids.add(value)
    return ids


def format_examples(values, limit=10):
    examples = sorted(values)[:limit]
    return ",".join(examples)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--feature-type", default="gene")
    parser.add_argument("--attribute", default="ID")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    count_ids = read_count_ids(args.counts)
    annotation_ids = read_annotation_ids(
        args.annotation, args.feature_type, args.attribute
    )
    overlap = count_ids & annotation_ids
    count_only = count_ids - annotation_ids
    annotation_only = annotation_ids - count_ids

    if not annotation_ids:
        raise SystemExit(
            "No annotation IDs found with "
            f"feature_type='{args.feature_type}' and attribute='{args.attribute}'."
        )

    count_fraction = len(overlap) / len(count_ids) if count_ids else 0
    annotation_fraction = len(overlap) / len(annotation_ids)

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["count_matrix_ids", len(count_ids)])
        writer.writerow(["annotation_ids", len(annotation_ids)])
        writer.writerow(["overlap_ids", len(overlap)])
        writer.writerow(["count_matrix_overlap_fraction", f"{count_fraction:.4f}"])
        writer.writerow(["annotation_overlap_fraction", f"{annotation_fraction:.4f}"])
        writer.writerow(["count_ids_without_annotation", len(count_only)])
        writer.writerow(["annotation_ids_without_counts", len(annotation_only)])
        writer.writerow(["count_only_examples", format_examples(count_only)])
        writer.writerow(["annotation_only_examples", format_examples(annotation_only)])

    print(
        "Count/annotation ID overlap: "
        f"{len(overlap)} of {len(count_ids)} count IDs overlap "
        f"{len(annotation_ids)} annotation IDs "
        f"(feature_type='{args.feature_type}', attribute='{args.attribute}').",
        file=sys.stderr,
    )

    if not overlap:
        raise SystemExit(
            "Count matrix IDs do not overlap annotation IDs. Check annotation "
            "feature_type/attribute settings and quantification output IDs."
        )


if __name__ == "__main__":
    main()
