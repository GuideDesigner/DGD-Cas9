#!/usr/bin/env python3
"""
sequence_utils.py — Shared DNA/FASTA sequence utilities for DGD
================================================================
Author : Vipin Menon, BIG Lab, Hanyang University
Date   : August 2021 (modernized 2026)

Provides:
  - reverse_complement()  — reverse-complement a DNA string
  - parse_fasta()         — lazy FASTA file parser
"""

import logging
from typing import Iterator, Optional, Tuple

logger = logging.getLogger(__name__)

# Precomputed DNA complement translation table (built once at import time)
_COMPLEMENT_TABLE = str.maketrans("ATCG", "TAGC")


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Args:
        sequence: A DNA string containing only A, T, C, G characters (uppercase).

    Returns:
        The reverse complement as an uppercase string.

    Example:
        >>> reverse_complement("ATCG")
        'CGAT'
    """
    return sequence[::-1].translate(_COMPLEMENT_TABLE)


def parse_fasta(filepath: str) -> Iterator[Tuple[str, str]]:
    """
    Parse a FASTA file and yield (sequence_id, sequence) pairs.

    Multi-line sequences are concatenated. Sequences are uppercased and
    whitespace-stripped. Empty lines are ignored.

    Args:
        filepath: Path to a FASTA-formatted file (.fa, .fasta, .fna).

    Yields:
        Tuples of (sequence_id, full_sequence).

    Raises:
        FileNotFoundError : If the file does not exist.
        ValueError        : If the file contains no valid FASTA records.

    Example:
        >>> for seq_id, sequence in parse_fasta("input.fa"):
        ...     print(seq_id, len(sequence))
    """
    import os
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"FASTA file not found: '{filepath}'")

    current_id: Optional[str] = None
    parts: list = []
    found_any = False

    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Flush previous record before starting a new one
                if current_id is not None:
                    yield current_id, "".join(parts).upper()
                    found_any = True
                # Take only the first word as the sequence ID
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

    # Flush the final record
    if current_id is not None:
        yield current_id, "".join(parts).upper()
        found_any = True

    if not found_any:
        raise ValueError(f"No valid FASTA records found in '{filepath}'")
