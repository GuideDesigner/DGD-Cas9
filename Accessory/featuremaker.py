#!/usr/bin/env python3
"""
Accessory/featuremaker.py — Standalone pipeline step 10
========================================================
Build per-region structural connectivity feature matrices and merge with
sequence features.
Usage: python featuremaker.py
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import build_features


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build structural feature matrix (step 10 of DGD pipeline)."
    )
    parser.add_argument(
        "--bp-input", default="spacer_scaffold_basepairs.csv",
        help="Spacer–scaffold base-pair CSV (default: spacer_scaffold_basepairs.csv)"
    )
    parser.add_argument(
        "--seq-input", default="Target_sequence_feature.csv",
        help="Target sequence features CSV (default: Target_sequence_feature.csv)"
    )
    parser.add_argument(
        "--annot-input", default="Structural_annotation.csv",
        help="Structural annotation CSV (default: Structural_annotation.csv)"
    )
    parser.add_argument(
        "--output", default="Feature_Data_Spacer_Scaffold.csv",
        help="Output combined feature CSV (default: Feature_Data_Spacer_Scaffold.csv)"
    )
    args = parser.parse_args()
    build_features(args.bp_input, args.seq_input, args.annot_input, args.output)


if __name__ == "__main__":
    main()
