#!/usr/bin/env python3
"""
Accessory/fastamaker.py — Standalone pipeline step 2
=====================================================
Append the sgRNA scaffold sequence to each guide candidate to prepare
sequences for RNAfold structure prediction.
Usage: python fastamaker.py [--input Structure_file.csv] [--output Structure_Connection.fa]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DGD import make_fasta_for_rnafold


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build scaffold-appended FASTA for RNAfold (step 2 of DGD pipeline)."
    )
    parser.add_argument(
        "--input", default="Structure_file.csv",
        help="Guide candidate CSV (default: Structure_file.csv)"
    )
    parser.add_argument(
        "--fasta-out", default="Structure_Connection.fa",
        help="Output FASTA file (default: Structure_Connection.fa)"
    )
    parser.add_argument(
        "--csv-out", default="Structure_Connection.csv",
        help="Output CSV file (default: Structure_Connection.csv)"
    )
    args = parser.parse_args()
    make_fasta_for_rnafold(args.input, args.fasta_out, args.csv_out)


if __name__ == "__main__":
    main()
