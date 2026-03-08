#!/usr/bin/env python3
"""
Accessory/score_deep.py — Standalone pipeline step 12
======================================================
Score guide candidates using the DGD ensemble of CNN models.
Usage: python score_deep.py [--models ./models] [--input Deep_learning_file.csv] [--output DGD.csv]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import score_guides


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Score guides with DGD ensemble (step 12 of DGD pipeline)."
    )
    parser.add_argument(
        "--models", default="./models",
        help="Directory containing trained .h5 model files (default: ./models)"
    )
    parser.add_argument(
        "--input", default="Deep_learning_file.csv",
        help="Deep learning input CSV (default: Deep_learning_file.csv)"
    )
    parser.add_argument(
        "--struct-input", default="Structure_file.csv",
        help="Guide structure CSV (default: Structure_file.csv)"
    )
    parser.add_argument(
        "--output", default="DGD.csv",
        help="Output scores CSV (default: DGD.csv)"
    )
    args = parser.parse_args()
    score_guides(args.models, args.input, args.struct_input, args.output)


if __name__ == "__main__":
    main()
