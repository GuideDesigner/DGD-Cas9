#!/usr/bin/env python3
"""
Accessory/serialconnection.py — Standalone pipeline step 9
===========================================================
Annotate each spacer–scaffold base pair with a structural region label
(R, TL, AR, LR, SL1, SL2, SL3, or NS).
Usage: python serialconnection.py [--input spacer_scaffold_feature.csv] [--output Structural_annotation.csv]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import annotate_structure_regions


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate structural regions (step 9 of DGD pipeline)."
    )
    parser.add_argument(
        "--input", default="spacer_scaffold_feature.csv",
        help="Connection frequency CSV (default: spacer_scaffold_feature.csv)"
    )
    parser.add_argument(
        "--output", default="Structural_annotation.csv",
        help="Output annotation CSV (default: Structural_annotation.csv)"
    )
    args = parser.parse_args()
    annotate_structure_regions(args.input, args.output)


if __name__ == "__main__":
    main()
