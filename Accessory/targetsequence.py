#!/usr/bin/env python3
"""
Accessory/targetsequence.py — Standalone pipeline step 3
=========================================================
Compute one-hot + sequence features for each guide candidate.
Usage: python targetsequence.py [--input Structure_file.csv] [--output Target_sequence_feature.csv]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import compute_target_features


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute target sequence features (step 3 of DGD pipeline)."
    )
    parser.add_argument(
        "--input", default="Structure_file.csv",
        help="Guide candidate CSV (default: Structure_file.csv)"
    )
    parser.add_argument(
        "--output", default="Target_sequence_feature.csv",
        help="Output feature CSV (default: Target_sequence_feature.csv)"
    )
    args = parser.parse_args()
    compute_target_features(args.input, args.output)


if __name__ == "__main__":
    main()
