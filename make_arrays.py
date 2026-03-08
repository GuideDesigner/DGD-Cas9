#!/usr/bin/env python3
"""
make_arrays.py — RNA structure array utilities for DGD
======================================================
Author : Vipin Menon, BIG Lab, Hanyang University
Date   : August 2021 (modernized 2026)

Provides helper functions for parsing and manipulating RNA dot-bracket
structures into partner-index arrays, dimer arrays, and connectivity arrays.
"""

from typing import List, Tuple


def return_partner_index(structure: str) -> List[int]:
    """
    Parse dot-bracket notation and return partner indices for each position.

    For each character in the dot-bracket string:
    - '.' → 0 (unpaired)
    - '(' → index of matching ')'
    - ')' → index of matching '('

    Args:
        structure: A dot-bracket string representing RNA secondary structure.

    Returns:
        List of integers where index i contains the partner position of base i,
        or 0 if the base is unpaired.
    """
    partner_index = []
    stack = []

    for i, char in enumerate(structure):
        if char == '.':
            partner_index.append(0)
        elif char == '(':
            stack.append(i)
            partner_index.append(0)  # placeholder
        elif char == ')':
            opening_index = stack.pop()
            partner_index[opening_index] = i + 1  # 1-indexed
            partner_index.append(opening_index + 1)  # 1-indexed

    return partner_index


def return_partner_base(partner_index: List[int], sequence: str) -> List[str]:
    """
    Return the base paired to each position, or '.' if unpaired.

    Args:
        partner_index: List of partner indices (from return_partner_index).
        sequence: RNA sequence string.

    Returns:
        List where index i contains the paired base or '.'.
    """
    partner_base = []
    for i in range(len(partner_index)):
        if partner_index[i] == 0:
            partner_base.append('.')
        else:
            partner_base.append(sequence[partner_index[i] - 1])
    return partner_base


def return_dimers_array(sequence: str, partner_index: List[int]) -> Tuple[List[str], List[str]]:
    """
    Extract consecutive base pairs (dimers) from sequence and partner index.

    For each position i in range(19), if both positions i and i+1 are paired,
    append the dimer and its complement.

    Args:
        sequence: RNA sequence string.
        partner_index: List of partner indices.

    Returns:
        Tuple of (dimers_list, anti_sequence_list).
    """
    dimers_array = []
    anti_seq = []

    for i in range(19):
        if partner_index[i] != 0 and partner_index[i + 1] != 0:
            dimers_array.append(sequence[i : i + 2])
            anti_seq.append(
                sequence[partner_index[i + 1] - 1]
                + sequence[partner_index[i] - 1]
            )

    return dimers_array, anti_seq


def adjust_partner_index(partner_index: List[int]) -> List[int]:
    """
    Adjust partner indices by zeroing isolated or short base pairs.

    Sets partner_index[i] = 0 if:
    - Both neighbors (i-1 and i+1) are unpaired, OR
    - partner_index[i] < 21

    Args:
        partner_index: List of partner indices.

    Returns:
        Modified partner_index list.
    """
    for i in range(1, 19):
        if (partner_index[i - 1] == 0 and partner_index[i + 1] == 0) or \
           partner_index[i] < 21:
            partner_index[i] = 0
    return partner_index


def return_masked_array(
    input_array: List[int],
    partner_index: List[int],
    size: int
) -> List[str]:
    """
    Create a masked array representation of partner indices.

    Args:
        input_array: Unused input array (for compatibility).
        partner_index: List of partner indices.
        size: Number of positions to process.

    Returns:
        List of strings: partner index values or ' ' for unpaired positions.
    """
    masked_array = []
    for i in range(size):
        if partner_index[i] == 0:
            masked_array.append(' ')
        else:
            masked_array.append(str(partner_index[i]))
    return masked_array


def return_connection_length(partner_index: List[int]) -> Tuple[List[int], int]:
    """
    Calculate lengths of consecutive paired-base runs.

    Iterates through partner_index, counting consecutive non-zero positions.

    Args:
        partner_index: List of partner indices.

    Returns:
        Tuple of (list of run lengths, total count of paired bases).
    """
    length_list = []
    total_count = 0
    i = 0

    while i < len(partner_index):
        if partner_index[i] != 0:
            length = 0
            while i < len(partner_index) and partner_index[i] != 0:
                length += 1
                total_count += 1
                i += 1
            length_list.append(length)
        else:
            i += 1

    return length_list, total_count


def convert_array_to_string(array: List[str]) -> str:
    """
    Join array elements into a single string.

    Args:
        array: List of strings.

    Returns:
        Concatenated string.
    """
    return ''.join(array)


def strip_string_into_array(string: str) -> List[str]:
    """
    Convert a stripped string into a list of characters.

    Args:
        string: Input string.

    Returns:
        List of characters.
    """
    return list(string.strip())


def remove_smaller_than_size(array: List[str], size: int) -> List[str]:
    """
    Filter array elements by minimum length.

    Args:
        array: List of strings to filter.
        size: Minimum length threshold.

    Returns:
        List of strings with length >= size.
    """
    return [element for element in array if len(element) >= size]
