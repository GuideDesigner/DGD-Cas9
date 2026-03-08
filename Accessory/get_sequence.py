#!/usr/bin/env python3
"""
Accessory/get_sequence.py — Standalone reverse complement and FASTA utilities
=============================================================================
Standalone access to shared DNA/FASTA utilities.
Usage: python get_sequence.py --sequence ATCGATCG
       python get_sequence.py --fasta input.fa
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from sequence_utils import reverse_complement, parse_fasta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reverse complement a sequence or parse a FASTA file."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sequence", "-s", help="DNA sequence to reverse-complement")
    group.add_argument("--fasta", "-f", help="FASTA file to parse")
    args = parser.parse_args()

    if args.sequence:
        print(reverse_complement(args.sequence.upper()))
    else:
        for seq_id, sequence in parse_fasta(args.fasta):
            print(f"{seq_id}\t{len(sequence)} nt")


if __name__ == "__main__":
    main()
