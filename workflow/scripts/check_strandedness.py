#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path


SKIP_ROWS = {
    "N_unmapped",
    "N_multimapping",
    "N_noFeature",
    "N_ambiguous",
}


def configured_strandedness(quantification_tool, star_value, featurecounts_value):
    if quantification_tool == "featurecounts":
        value = str(featurecounts_value)
        return {"0": "unstranded", "1": "forward", "2": "reverse"}.get(
            value, value
        )
    return str(star_value).lower().replace("-", "_")


def infer_star_strandedness(path, min_fraction, min_ratio):
    unstranded = 0
    forward = 0
    reverse = 0

    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] in SKIP_ROWS:
                continue
            unstranded += int(fields[1])
            forward += int(fields[2])
            reverse += int(fields[3])

    stranded_total = forward + reverse
    if stranded_total == 0:
        return "unknown", unstranded, forward, reverse, 0.0, 0.0

    dominant = max(forward, reverse)
    minor = min(forward, reverse)
    fraction = dominant / stranded_total
    ratio = float("inf") if minor == 0 else dominant / minor

    if fraction >= min_fraction and ratio >= min_ratio:
        inferred = "forward" if forward > reverse else "reverse"
    else:
        inferred = "unstranded"

    return inferred, unstranded, forward, reverse, fraction, ratio


def main():
    parser = argparse.ArgumentParser(
        description="Infer library strandedness from STAR ReadsPerGene output."
    )
    parser.add_argument("--star-counts", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--quantification-tool", required=True)
    parser.add_argument("--star-strandedness", required=True)
    parser.add_argument("--featurecounts-strandedness", required=True)
    parser.add_argument("--mode", choices=["warn", "fail"], default="warn")
    parser.add_argument("--min-fraction", type=float, default=0.75)
    parser.add_argument("--min-ratio", type=float, default=2.0)
    args = parser.parse_args()

    quantification_tool = args.quantification_tool.lower()
    expected = configured_strandedness(
        quantification_tool,
        args.star_strandedness,
        args.featurecounts_strandedness,
    )
    inferred, unstranded, forward, reverse, fraction, ratio = infer_star_strandedness(
        args.star_counts,
        args.min_fraction,
        args.min_ratio,
    )

    status = "pass"
    message = ""
    if inferred != "unknown" and expected != inferred:
        status = "mismatch"
        message = (
            f"Configured strandedness is '{expected}', but STAR gene counts "
            f"suggest '{inferred}' for sample '{args.sample}'."
        )
    elif inferred == "unknown":
        status = "unknown"
        message = f"Unable to infer strandedness for sample '{args.sample}'."

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample",
                "configured_strandedness",
                "inferred_strandedness",
                "status",
                "star_unstranded_counts",
                "star_forward_counts",
                "star_reverse_counts",
                "dominant_strand_fraction",
                "dominant_to_minor_ratio",
                "message",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample": args.sample,
                "configured_strandedness": expected,
                "inferred_strandedness": inferred,
                "status": status,
                "star_unstranded_counts": unstranded,
                "star_forward_counts": forward,
                "star_reverse_counts": reverse,
                "dominant_strand_fraction": f"{fraction:.4f}",
                "dominant_to_minor_ratio": (
                    "inf" if ratio == float("inf") else f"{ratio:.4f}"
                ),
                "message": message,
            }
        )

    if status == "mismatch":
        print(message, file=sys.stderr)
        if args.mode == "fail":
            return 1
    if status == "unknown":
        print(message, file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
