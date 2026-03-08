#!/usr/bin/env python3
"""
Accessory/finalfeatures.py — Standalone pipeline step 11
=========================================================
Compute final spacer–scaffold structural features (GC ratio, Gibbs energy,
monomer and dimer counts) and merge with feature data.
Usage: python finalfeatures.py
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import compute_final_features


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute final structural features (step 11 of DGD pipeline)."
    )
    parser.add_argument(
        "--rnafold-out", default="Structure_Connection.out",
        help="RNAfold dot-bracket output (default: Structure_Connection.out)"
    )
    parser.add_argument(
        "--feature-input", default="Feature_Data_Spacer_Scaffold.csv",
        help="Combined feature CSV (default: Feature_Data_Spacer_Scaffold.csv)"
    )
    parser.add_argument(
        "--output", default="Deep_learning_file.csv",
        help="Output deep learning input CSV (default: Deep_learning_file.csv)"
    )
    args = parser.parse_args()
    compute_final_features(args.rnafold_out, args.feature_input, args.output)


if __name__ == "__main__":
    main()
