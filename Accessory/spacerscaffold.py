#!/usr/bin/env python3
"""
Accessory/spacerscaffold.py — Standalone pipeline step 7
=========================================================
Extract spacer–scaffold base-pair columns from the full connection matrix.
Usage: python spacerscaffold.py [--input Structure_basepairs.csv] [--output spacer_scaffold_basepairs.csv]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import extract_spacer_scaffold_pairs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract spacer–scaffold base pairs (step 7 of DGD pipeline)."
    )
    parser.add_argument(
        "--input", default="Structure_basepairs.csv",
        help="Full connection matrix CSV (default: Structure_basepairs.csv)"
    )
    parser.add_argument(
        "--output", default="spacer_scaffold_basepairs.csv",
        help="Output spacer–scaffold CSV (default: spacer_scaffold_basepairs.csv)"
    )
    args = parser.parse_args()
    extract_spacer_scaffold_pairs(args.input, args.output)


if __name__ == "__main__":
    main()
