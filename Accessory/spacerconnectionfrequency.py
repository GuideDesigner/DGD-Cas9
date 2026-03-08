#!/usr/bin/env python3
"""
Accessory/spacerconnectionfrequency.py — Standalone pipeline step 8
====================================================================
Build connection-frequency position metadata from spacer–scaffold base pairs.
Usage: python spacerconnectionfrequency.py [--input spacer_scaffold_basepairs.csv] [--output spacer_scaffold_feature.csv]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import compute_connection_frequency


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build connection frequency table (step 8 of DGD pipeline)."
    )
    parser.add_argument(
        "--input", default="spacer_scaffold_basepairs.csv",
        help="Spacer–scaffold base-pair CSV (default: spacer_scaffold_basepairs.csv)"
    )
    parser.add_argument(
        "--output", default="spacer_scaffold_feature.csv",
        help="Output feature CSV (default: spacer_scaffold_feature.csv)"
    )
    args = parser.parse_args()
    compute_connection_frequency(args.input, args.output)


if __name__ == "__main__":
    main()
