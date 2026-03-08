#!/usr/bin/env python3
"""
Accessory/Connectstructure.py — Standalone pipeline step 5
===========================================================
Parse b2ct-formatted RNAfold output into a tabular connection matrix.
Usage: python Connectstructure.py [--input Structure_Connection.outs] [--output Structure_out.txt]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import parse_rnafold_output


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse RNAfold b2ct output into connection matrix (step 5 of DGD pipeline)."
    )
    parser.add_argument(
        "--input", default="Structure_Connection.outs",
        help="b2ct output file (default: Structure_Connection.outs)"
    )
    parser.add_argument(
        "--output", default="Structure_out.txt",
        help="Output connection matrix (default: Structure_out.txt)"
    )
    args = parser.parse_args()
    parse_rnafold_output(args.input, args.output)


if __name__ == "__main__":
    main()
