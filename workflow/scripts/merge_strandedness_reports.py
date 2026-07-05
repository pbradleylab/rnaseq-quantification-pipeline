#!/usr/bin/env python3

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-sample strandedness check reports."
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("reports", nargs="+")
    args = parser.parse_args()

    header_written = False
    with open(args.output, "w", newline="") as output_handle:
        writer = None
        for report in args.reports:
            with open(report, newline="") as input_handle:
                reader = csv.DictReader(input_handle, delimiter="\t")
                if writer is None:
                    writer = csv.DictWriter(
                        output_handle,
                        fieldnames=reader.fieldnames,
                        delimiter="\t",
                    )
                if not header_written:
                    writer.writeheader()
                    header_written = True
                for row in reader:
                    writer.writerow(row)


if __name__ == "__main__":
    main()
