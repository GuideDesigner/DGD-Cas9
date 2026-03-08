#!/usr/bin/env python3
"""
DGD — Deep Guide Designer
==========================
Copyright : Vipin Menon, Jang-il Sohn, Seokju Park & BIG Lab, Hanyang University (HYU)
Author    : Vipin Menon  <a.vipin.menon@gmail.com>
Date      : 21 August 2021 (modernized 2026)

Description
-----------
DGD predicts CRISPR-Cas9 sgRNA on-target activity by integrating three
information sources into a CNN-based deep learning model:

  1. Target sequence features   (nucleotide composition, free energy, GC content)
  2. Spacer–scaffold base-pairs (structural connectivity features)
  3. Sequence context           (one-hot encoded 30-bp window)

Input  : FASTA file (.fa) with one or more sequences (100–10 000 nt each)
Output : DGD.csv  —  ID, Start, End, Strand, Sequence, DGD score

Usage
-----
  python DGD.py input.fa
  python DGD.py input.fa --output results.csv --models ./models --verbose
"""

# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import argparse
import itertools
import logging
import os
import subprocess
import sys
from collections import OrderedDict, defaultdict
from itertools import chain
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Third-party imports (with clear error messages if missing)
# ---------------------------------------------------------------------------
import numpy as np
import pandas as pd

try:
    import RNA
except ImportError as exc:
    raise ImportError(
        "ViennaRNA Python bindings are required.\n"
        "Install with:  conda install -c bioconda viennarna"
    ) from exc

try:
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError as exc:
    raise ImportError(
        "Biopython is required.\n"
        "Install with:  pip install biopython  or  conda install biopython"
    ) from exc

try:
    import tensorflow as tf
    from tensorflow.keras import models as keras_models
except ImportError as exc:
    raise ImportError(
        "TensorFlow 2.x is required.\n"
        "Install with:  pip install tensorflow==2.2.0"
    ) from exc

import make_arrays
import stacking_model
from sequence_utils import parse_fasta, reverse_complement

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pipeline constants
# ---------------------------------------------------------------------------

#: sgRNA scaffold sequence appended during FASTA generation for RNAfold
SCAFFOLD_SEQ: str = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAA"
    "AGTGGCACCGAGTCGGTGCTTTTTT"
)

#: Nucleotide alphabet
NUCLEOTIDES: List[str] = ["A", "T", "G", "C"]

#: All dinucleotides in alphabetical order
DINUCLEOTIDES: List[str] = [a + b for b in NUCLEOTIDES for a in NUCLEOTIDES]

#: Structural region labels and their scaffold position ranges
STRUCTURE_REGIONS: Dict[str, Tuple[int, int]] = {
    "R":   (21, 32),
    "TL":  (33, 36),
    "AR":  (37, 49),
    "LR":  (63, 67),
    "SL1": (54, 58),
    "SL2": (73, 76),
    "SL3": (88, 90),
    "NS":  (0,  0),   # "Not Specified" — filled by default after labelling others
}


class PipelineFiles:
    """
    Central registry for all intermediate file paths used by the DGD pipeline.

    Keeping paths in one place makes it easy to redirect I/O, run the pipeline
    in a custom working directory, or swap in different file names without
    hunting through the codebase.
    """
    STRUCTURE_FILE:     str = "Structure_file.csv"
    TARGET_FEATURES:    str = "Target_sequence_feature.csv"
    FASTA_OUT:          str = "Structure_Connection.fa"
    SEQUENCE_CSV:       str = "Structure_Connection.csv"
    RNAFOLD_OUT:        str = "Structure_Connection.out"
    RNAFOLD_OUTS:       str = "Structure_Connection.outs"
    STRUCTURE_OUT_TMP:  str = "Structure_Cas9_out.txt"
    STRUCTURE_OUT:      str = "Structure_out.txt"
    BASEPAIRS_CSV:      str = "Structure_basepairs.csv"
    SPACER_CSV:         str = "spacer_scaffold_basepairs.csv"
    SPACER_FEATURE_CSV: str = "spacer_scaffold_feature.csv"
    STRUCTURAL_ANNOT:   str = "Structural_annotation.csv"
    FEATURE_DATA:       str = "Feature_Data_Spacer_Scaffold.csv"
    DEEP_LEARNING:      str = "Deep_learning_file.csv"
    OUTPUT:             str = "DGD.csv"


# ---------------------------------------------------------------------------
# Step 1 — Scan FASTA and write guide candidates to Structure_file.csv
# ---------------------------------------------------------------------------

def scan_guides(fasta_path: str) -> None:
    """
    Scan a FASTA file for all Cas9 guide RNA candidates (NGG PAM) on both
    strands and write them to ``Structure_file.csv``.

    Each 30-bp window centred on a GG (or CC on the reverse strand) PAM is
    scored for position validity and saved as a candidate guide.

    Args:
        fasta_path: Path to the input FASTA file (sequences 100–10 000 nt).

    Raises:
        FileNotFoundError : If the FASTA file does not exist.
        ValueError        : If no valid sequences are found.
    """
    forward_guides: Dict[str, List] = {}
    reverse_guides: Dict[str, List] = {}

    for seq_id, sequence in parse_fasta(fasta_path):
        seq_len = len(sequence)
        if not (100 <= seq_len <= 10_000):
            logger.warning(
                "Sequence '%s' length %d is outside the 100–10 000 nt range. Skipping.",
                seq_id, seq_len,
            )
            continue

        # Forward strand: PAM = NGG (look for GG, take 25 bp upstream + 5 bp)
        search_pos = 0
        while True:
            pos = sequence.find("GG", search_pos)
            if pos == -1:
                break
            if 25 < pos < seq_len - 5:
                start, end = pos - 25, pos + 5
                guide_id = f"{seq_id}:{start}:{end}"
                forward_guides[guide_id] = [start, end, "+", sequence[start:end]]
            search_pos = pos + 1

        # Reverse strand: PAM = NCC on forward ≡ NGG on reverse
        search_pos = 0
        while True:
            pos = sequence.find("CC", search_pos)
            if pos == -1:
                break
            if 3 < pos < seq_len - 27:
                start, end = pos - 3, pos + 27
                rc_window = reverse_complement(sequence[start:end])
                guide_id = f"{seq_id}:{start}:{end}"
                reverse_guides[guide_id] = [start, end, "-", rc_window]
            search_pos = pos + 1

    all_guides = dict(chain(forward_guides.items(), reverse_guides.items()))
    sorted_guides = OrderedDict(
        sorted(all_guides.items(), key=lambda item: item[1][3], reverse=True)
    )

    with open(PipelineFiles.STRUCTURE_FILE, "w", encoding="utf-8") as out:
        out.write("ID,Start,End,Strand,Sequence\n")
        for guide_id, (start, end, strand, seq) in sorted_guides.items():
            out.write(f"{guide_id},{start},{end},{strand},{seq}\n")

    logger.info(
        "Guide scan complete: %d candidates written to '%s'",
        len(sorted_guides), PipelineFiles.STRUCTURE_FILE,
    )


# ---------------------------------------------------------------------------
# Step 2 — Generate FASTA for RNAfold (spacer + scaffold)
# ---------------------------------------------------------------------------

def make_fasta_for_rnafold() -> None:
    """
    Append the Cas9 scaffold sequence to each guide's spacer region and write
    a FASTA file for downstream RNAfold structure prediction.

    Reads  : ``Structure_file.csv``
    Writes : ``Structure_Connection.fa``, ``Structure_Connection.csv``
    """
    struct_df = pd.read_csv(PipelineFiles.STRUCTURE_FILE)

    with (
        open(PipelineFiles.FASTA_OUT, "w", encoding="utf-8") as fa_out,
        open(PipelineFiles.SEQUENCE_CSV, "w", encoding="utf-8") as csv_out,
    ):
        csv_out.write("ID,Sequence\n")
        for _, row in struct_df.iterrows():
            guide_id = row["ID"]
            spacer   = row["Sequence"][4:24]          # 20-bp spacer (skip PAM)
            full_seq = spacer + SCAFFOLD_SEQ           # spacer + scaffold

            fa_out.write(f">{guide_id}\n{full_seq}\n")
            csv_out.write(f"{guide_id},{full_seq}\n")

    logger.info("RNAfold FASTA written to '%s'", PipelineFiles.FASTA_OUT)


# ---------------------------------------------------------------------------
# Step 3 — Compute target sequence features
# ---------------------------------------------------------------------------

def _build_position_features(sequence: str) -> Tuple[List[int], List[int]]:
    """
    Build per-position one-hot encoding vectors for mono- and dinucleotides.

    Args:
        sequence: The full 30-bp guide sequence (uppercase).

    Returns:
        Tuple of (mononucleotide_vector, dinucleotide_vector),
        each as a flat list of 0/1 integers.
    """
    mono = [
        1 if nuc == base else 0
        for nuc in NUCLEOTIDES
        for base in sequence
    ]
    di = [
        1 if dinuc == sequence[i:i+2] else 0
        for dinuc in DINUCLEOTIDES
        for i in range(len(sequence) - 1)
    ]
    return mono, di


def _sequence_entropy(region: str) -> float:
    """
    Compute Shannon entropy of nucleotide composition in *region*.

    Args:
        region: A nucleotide string (uppercase).

    Returns:
        Rounded Shannon entropy (bits).
    """
    length = len(region)
    entropy = 0.0
    for nuc in NUCLEOTIDES:
        freq = region.count(nuc) / float(length)
        if freq > 0:
            entropy -= freq * np.log2(freq)
    return round(entropy, 1)


def compute_target_features() -> None:
    """
    Extract sequence, thermodynamic, and compositional features for every guide
    candidate and write them to ``Target_sequence_feature.csv``.

    Reads  : ``Structure_file.csv``
    Writes : ``Target_sequence_feature.csv``

    Features computed per guide:
      - Per-position mono/dinucleotide one-hot encoding
      - Shannon entropy of the 20-bp target region
      - RNA free energy (ViennaRNA)
      - GC content percentage, high/low binary flags
      - Melting temperature (nearest-neighbour method)
      - Global mono/dinucleotide counts
    """
    # Build column header
    mono_cols  = [f"{n}{i}" for n in NUCLEOTIDES for i in range(1, 31)]
    di_cols    = [f"{d}{i}" for d in DINUCLEOTIDES for i in range(1, 30)]
    extra_cols = [
        "Entropy", "Energy", "GCcount", "Gchigh", "GClow",
        "MeltingTemperature",
        "A", "T", "G", "C",
        "AA", "AT", "AG", "AC", "CA", "CG", "CC", "CT",
        "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT",
    ]
    header = ["ID"] + mono_cols + di_cols + extra_cols

    struct_df = pd.read_csv(PipelineFiles.STRUCTURE_FILE)
    rows: List[Dict] = []

    for _, row in struct_df.iterrows():
        guide_id  = row["ID"]
        full_seq  = row["Sequence"].upper()
        target    = full_seq[4:24]     # 20-bp target region (positions 4–23)

        mono_vec, di_vec = _build_position_features(full_seq)

        # Thermodynamic features
        entropy      = _sequence_entropy(target)
        free_energy  = round(RNA.fold(target)[-1], 0)
        gc_count     = target.count("G") + target.count("C")
        gc_content   = round(gc_count / float(len(target)) * 100, 0)
        gc_high      = 1 if gc_count >= 10 else 0
        gc_low       = 1 if gc_count < 10  else 0
        melt_temp    = mt.Tm_NN(target)

        # Global nucleotide / dinucleotide counts
        global_counts = {nuc: full_seq.count(nuc) for nuc in NUCLEOTIDES}
        global_di     = {di: full_seq.count(di) for di in DINUCLEOTIDES}

        record: Dict = {"ID": guide_id}
        for i, val in enumerate(mono_vec):
            record[header[1 + i]] = val
        offset = 1 + len(mono_vec)
        for i, val in enumerate(di_vec):
            record[header[offset + i]] = val

        record.update({
            "Entropy":           entropy,
            "Energy":            free_energy,
            "GCcount":           gc_content,
            "Gchigh":            gc_high,
            "GClow":             gc_low,
            "MeltingTemperature": melt_temp,
            **{n: global_counts[n] for n in NUCLEOTIDES},
            **{d: global_di[d]     for d in DINUCLEOTIDES},
        })
        rows.append(record)

    pd.DataFrame(rows, columns=header).to_csv(
        PipelineFiles.TARGET_FEATURES, index=False
    )
    logger.info(
        "Target features written to '%s'  (%d guides)",
        PipelineFiles.TARGET_FEATURES, len(rows),
    )


# ---------------------------------------------------------------------------
# Step 4 — Parse RNAfold dot-bracket output and build connection matrix
# ---------------------------------------------------------------------------

def parse_rnafold_output() -> None:
    """
    Parse RNAfold dot-bracket output and build a tabular connection matrix.

    Reads  : ``Structure_Connection.outs``  (b2ct format)
    Writes : ``Structure_out.txt``

    Each row of the output represents one guide and contains 102 connection
    columns (one per position in spacer + scaffold).
    """
    header_cols = ["ID"] + [f"Pos{i}" for i in range(1, 103)]

    guide_connections: Dict[str, List[str]] = {}
    current_id: Optional[str] = None

    with open(PipelineFiles.RNAFOLD_OUTS, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = list(filter(None, line.strip().split(" ")))
            if len(parts) == 5:
                current_id = parts[4]
                guide_connections[current_id] = []
            elif current_id is not None and len(parts) >= 5:
                guide_connections[current_id].append(parts[4])

    df = pd.DataFrame.from_dict(guide_connections, orient="index")
    df.index.name = "ID"
    df.reset_index(inplace=True)
    df.columns = header_cols
    df.to_csv(PipelineFiles.STRUCTURE_OUT, sep="\t", index=False)

    logger.info("Connection matrix written to '%s'", PipelineFiles.STRUCTURE_OUT)


# ---------------------------------------------------------------------------
# Step 5 — Extract spacer–scaffold base-pair columns
# ---------------------------------------------------------------------------

def extract_spacer_scaffold_pairs() -> None:
    """
    Extract spacer (positions 1–20) to scaffold (positions 21–102) base-pair
    columns from the connection matrix.

    Reads  : ``Structure_basepairs.csv``
    Writes : ``spacer_scaffold_basepairs.csv``
    """
    conn_df = pd.read_csv(PipelineFiles.BASEPAIRS_CSV)

    pair_cols = [
        f"Connection_Pos{spacer}_Pos{scaffold}"
        for spacer   in range(1, 21)
        for scaffold in range(21, 103)
    ]
    pair_cols.append("ID")

    conn_df[pair_cols].to_csv(PipelineFiles.SPACER_CSV, index=False)
    logger.info(
        "Spacer–scaffold pairs written to '%s'", PipelineFiles.SPACER_CSV
    )


# ---------------------------------------------------------------------------
# Step 6 — Compute spacer–scaffold connection frequency
# ---------------------------------------------------------------------------

def compute_connection_frequency() -> None:
    """
    Pivot the spacer–scaffold connection table to a long format and extract
    position labels.

    Reads  : ``spacer_scaffold_basepairs.csv``
    Writes : ``spacer_scaffold_feature.csv``
    """
    conn_df = pd.read_csv(PipelineFiles.SPACER_CSV)
    long_df = conn_df.set_index("ID").T.reset_index()

    feat_df = pd.DataFrame(long_df.iloc[:, 0])
    feat_df.columns = ["nucleotide"]
    feat_df[["_label", "Pos_A", "Pos_B"]] = feat_df["nucleotide"].str.split(
        "_", n=2, expand=True
    ).iloc[:, [0, 1, 2]]
    feat_df["Pos_A"] = feat_df["Pos_A"].str.extract(r"(\d+)", expand=False).astype(int)
    feat_df["Pos_B"] = feat_df["Pos_B"].str.extract(r"(\d+)", expand=False).astype(int)
    feat_df = feat_df[["nucleotide", "Pos_A", "Pos_B"]]

    feat_df.to_csv(PipelineFiles.SPACER_FEATURE_CSV, index=False)
    logger.info(
        "Connection frequency written to '%s'", PipelineFiles.SPACER_FEATURE_CSV
    )


# ---------------------------------------------------------------------------
# Step 7 — Annotate structural regions
# ---------------------------------------------------------------------------

def annotate_structure_regions() -> None:
    """
    Label each spacer–scaffold connection by the scaffold structural region it
    maps to (R, TL, AR, LR, SL1, SL2, SL3, or NS for none of the above).

    Region definitions (scaffold positions):
      R   : 21–32    (repeat)
      TL  : 33–36    (tetra loop)
      AR  : 37–49    (anti-repeat)
      LR  : 63–67    (linking region)
      SL1 : 54–58    (stem loop 1)
      SL2 : 73–76    (stem loop 2)
      SL3 : 88–90    (stem loop 3)
      NS  : all others (not specified)

    Reads  : ``spacer_scaffold_feature.csv``
    Writes : ``Structural_annotation.csv``
    """
    df = pd.read_csv(PipelineFiles.SPACER_FEATURE_CSV)

    region_ranges = {
        "R":   (21, 32),
        "TL":  (33, 36),
        "AR":  (37, 49),
        "LR":  (63, 67),
        "SL1": (54, 58),
        "SL2": (73, 76),
        "SL3": (88, 90),
    }

    for label, (lo, hi) in region_ranges.items():
        mask = (df["Pos_B"] >= lo) & (df["Pos_B"] <= hi)
        df.loc[mask, "Structure"] = label

    df["Structure"].fillna("NS", inplace=True)
    df.to_csv(PipelineFiles.STRUCTURAL_ANNOT, index=False)
    logger.info("Structural annotation written to '%s'", PipelineFiles.STRUCTURAL_ANNOT)


# ---------------------------------------------------------------------------
# Step 8 — Build per-guide structural connectivity features
# ---------------------------------------------------------------------------

def _build_structure_df(
    region_df: pd.DataFrame,
    connection_bp: pd.DataFrame,
    prefix: str,
) -> pd.DataFrame:
    """
    Build a DataFrame of binary connectivity features for one structural region.

    For each unique spacer position in *region_df*, determine whether each
    guide has any base-pair connection to that position's scaffold partners.
    Sum all per-position flags to produce a ``Total_<prefix>`` column.

    This function eliminates the 8 near-identical copy-pasted blocks that
    existed in the original ``featuremaker()`` function.

    Args:
        region_df:     Rows from the structural annotation for one region label.
        connection_bp: Full spacer–scaffold base-pair DataFrame.
        prefix:        Region label used to name columns (e.g., "AR", "SL1").

    Returns:
        A DataFrame with one row per guide ID and one binary column per spacer
        position, plus a ``Total_<prefix>`` summary column.
    """
    positions = sorted(
        set(region_df["Pos_A"].tolist()),
        reverse=True,
    )
    nested: Dict[str, Dict[str, int]] = {}

    for pos in positions:
        col_name  = f"{prefix}{pos}"
        # Nucleotides that connect from this spacer position
        nuc_list  = list(set(
            region_df[region_df["Pos_A"] == pos]["nucleotide"].tolist()
        ))
        if not nuc_list:
            continue

        conn_info = connection_bp[nuc_list].copy()
        conn_info.insert(0, "ID", connection_bp["ID"])
        id_dict   = conn_info.set_index("ID").T.to_dict("list")

        pos_flags: Dict[str, int] = {}
        for guide_id, values in id_dict.items():
            pos_flags[guide_id] = 1 if sum(values) > 0 else 0

        nested[col_name] = pos_flags

    result_df = pd.DataFrame(nested)
    result_df[f"Total_{prefix}"] = result_df.sum(axis=1)
    result_df.index.name = "ID"
    result_df.reset_index(inplace=True)
    return result_df


def build_features() -> None:
    """
    Combine target sequence features with spacer–scaffold structural connectivity
    features into a single feature matrix for deep learning.

    Reads  : ``Structure_basepairs.csv``, ``Target_sequence_feature.csv``,
             ``Structural_annotation.csv``
    Writes : ``Feature_Data_Spacer_Scaffold.csv``
    """
    conn_bp      = pd.read_csv(PipelineFiles.BASEPAIRS_CSV)
    sequence_df  = pd.read_csv(PipelineFiles.TARGET_FEATURES)
    annot_df     = pd.read_csv(PipelineFiles.STRUCTURAL_ANNOT)

    # Build one connectivity DataFrame per structural region
    region_labels = ["AR", "NS", "R", "SL1", "LR", "SL2", "SL3", "TL"]
    region_dfs: Dict[str, pd.DataFrame] = {}

    for label in region_labels:
        subset = annot_df[annot_df["Structure"] == label][
            ["nucleotide", "Pos_A", "Pos_B", "Structure"]
        ]
        region_dfs[label] = _build_structure_df(subset, conn_bp, label)

    # Aggregate total connections across all regions
    total_conn = pd.DataFrame(
        pd.concat(
            [region_dfs[lbl][[f"Total_{lbl}"]] for lbl in region_labels],
            axis=1,
        ).sum(axis=1),
        columns=["Total_Connection"],
    )

    # Build the full connectivity table
    ar_df = region_dfs["AR"]
    other_frames = [
        region_dfs[lbl].drop(columns=["ID"])
        for lbl in ["R", "LR", "SL2", "TL", "NS", "SL1", "SL3"]
    ]
    connectivity = pd.concat([ar_df] + other_frames + [total_conn], axis=1)

    # Merge with sequence features
    merged = sequence_df.merge(connectivity, on="ID", how="inner")
    merged.to_csv(PipelineFiles.FEATURE_DATA, index=False)
    logger.info(
        "Combined features written to '%s'  (%d guides, %d features)",
        PipelineFiles.FEATURE_DATA, len(merged), len(merged.columns),
    )


# ---------------------------------------------------------------------------
# Step 9 — Compute spacer–scaffold monomer / dimer / energy features
# ---------------------------------------------------------------------------

def read_rnafold_dot_bracket(filepath: str) -> Dict[str, List[str]]:
    """
    Parse an RNAfold dot-bracket file (3 lines per entry: ID, sequence, structure).

    Args:
        filepath: Path to the ``.out`` file produced by RNAfold.

    Returns:
        Dict mapping guide ID → [sequence, dot-bracket structure, free energy string].
    """
    records: Dict[str, List[str]] = {}

    with open(filepath, "r", encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh]

    for i in range(0, len(lines) - 2, 3):
        guide_id  = lines[i].lstrip(">").strip()
        sequence  = lines[i + 1]
        struct_ln = lines[i + 2]
        structure = struct_ln[:struct_ln.rfind("(")]
        energy    = struct_ln[struct_ln.rfind("("):]
        records[guide_id] = [sequence, structure, energy]

    return records


def _count_monomers(sequence: str, partner_index: List[int]) -> Dict[str, int]:
    """
    Count monomers in paired positions of the spacer–scaffold junction.

    Args:
        sequence:      Spacer + scaffold nucleotide sequence.
        partner_index: Index list where partner_index[i] gives the pairing partner.

    Returns:
        Dict of nucleotide → count for paired bases.
    """
    counts: Dict[str, int] = {"A": 0, "C": 0, "G": 0, "U": 0}
    for i in range(20):
        if partner_index[i] != 0:
            counts[sequence[i]]                  += 1
            counts[sequence[partner_index[i]]]   += 1
    return counts


def _count_dimers(dimers: List[str]) -> Dict[str, int]:
    """
    Count di-nucleotide occurrences in a paired-region dimers list.

    Args:
        dimers: List of 2-character nucleotide strings.

    Returns:
        Dict of dinucleotide → count.
    """
    rna_dinucs = [a + b for a in "ACGU" for b in "ACGU"]
    counts: Dict[str, int] = {d: 0 for d in rna_dinucs}
    for dimer in dimers:
        if dimer in counts:
            counts[dimer] += 1
    return counts


def _count_gc_pairs(sequence: str, partner_index: List[int]) -> int:
    """
    Count G–C base pairs in the spacer region (positions 0–19).

    Args:
        sequence:      Full sequence string.
        partner_index: Pairing partner indices.

    Returns:
        Number of G–C (and C–G) Watson-Crick base pairs.
    """
    count = 0
    for i in range(20):
        if partner_index[i] != 0:
            pair = sequence[i] + sequence[partner_index[i]]
            if pair in ("GC", "CG"):
                count += 1
    return count


def _calculate_stacking_energy(dimers: List[str]) -> float:
    """
    Estimate free energy of a paired region from its dinucleotide stacking model.

    Args:
        dimers: Dinucleotide pair list from ``make_arrays.return_dimers_array``.

    Returns:
        Cumulative stacking free energy (kcal/mol).
    """
    stacking = stacking_model.return_stacking_model()
    energy   = 0.0
    for i in range(0, len(dimers), 2):
        energy += stacking[dimers[i]][dimers[i + 1][::-1]]
    return energy


def compute_final_features() -> None:
    """
    Compute spacer–scaffold monomer/dimer counts, GC ratio, and stacking energy
    for each guide from the RNAfold structural output, then merge with the
    combined feature matrix.

    Reads  : ``Structure_Connection.out``, ``Feature_Data_Spacer_Scaffold.csv``
    Writes : ``Deep_learning_file.csv``
    """
    structure_data = read_rnafold_dot_bracket(PipelineFiles.RNAFOLD_OUT)
    feature_df     = pd.read_csv(PipelineFiles.FEATURE_DATA)

    records_gc_energy: Dict[str, List] = {}
    records_counts:    Dict[str, Dict] = {}

    for guide_id, (sequence, structure, _energy) in structure_data.items():
        partner_idx   = make_arrays.return_partner_index(structure)
        partner_idx   = make_arrays.adjust_partner_index(partner_idx)
        partner_base  = make_arrays.return_partner_base(partner_idx, sequence)
        _conn_len, count = make_arrays.return_connection_length(partner_idx)
        dimers        = make_arrays.return_dimers_array(sequence, partner_idx)

        gc_ratio       = _count_gc_pairs(sequence, partner_idx) / count if count else 0.0
        monomer_counts = _count_monomers(sequence, partner_idx)
        dimer_counts   = _count_dimers(dimers)
        stacking_e     = _calculate_stacking_energy(dimers)

        records_gc_energy[guide_id] = [gc_ratio, stacking_e]
        records_counts[guide_id]    = {**monomer_counts, **dimer_counts}

    # Build DataFrames
    df_counts = pd.DataFrame(records_counts).T
    df_counts.index.name = "ID"
    df_counts.reset_index(inplace=True)

    df_energy = pd.DataFrame.from_dict(records_gc_energy, orient="index",
                                       columns=["GC_ratio", "Gibbs_Energy"])
    df_energy.index.name = "ID"
    df_energy.reset_index(inplace=True)

    scaffold_features = pd.merge(df_counts, df_energy, on="ID", how="inner")
    scaffold_features.columns = [
        "ID",
        "Spacer_Scaffold_A", "Spacer_Scaffold_C",
        "Spacer_Scaffold_G", "Spacer_Scaffold_U",
        "Spacer_Scaffold_AA", "Spacer_Scaffold_AC",
        "Spacer_Scaffold_AG", "Spacer_Scaffold_AU",
        "Spacer_Scaffold_CA", "Spacer_Scaffold_CC",
        "Spacer_Scaffold_CG", "Spacer_Scaffold_CU",
        "Spacer_Scaffold_GA", "Spacer_Scaffold_GC",
        "Spacer_Scaffold_GG", "Spacer_Scaffold_GU",
        "Spacer_Scaffold_UA", "Spacer_Scaffold_UC",
        "Spacer_Scaffold_UG", "Spacer_Scaffold_UU",
        "Spacer_Scaffold_GC_ratio", "Spacer_Scaffold_Gibbs_Energy",
    ]

    final_df = feature_df.merge(scaffold_features, on="ID", how="inner")
    final_df.to_csv(PipelineFiles.DEEP_LEARNING, index=False)
    logger.info(
        "Deep learning features written to '%s'  (%d guides, %d features)",
        PipelineFiles.DEEP_LEARNING, len(final_df), len(final_df.columns),
    )


# ---------------------------------------------------------------------------
# Step 10 — Prepare model inputs
# ---------------------------------------------------------------------------

def prepare_model_inputs(
    feature_csv: str,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Load the deep learning feature CSV and convert it to model-ready arrays.

    Returns a one-hot sequence tensor (shape: N×30×4×1), an auxiliary feature
    matrix, and the list of guide IDs.

    Args:
        feature_csv: Path to the deep learning feature CSV file.

    Returns:
        Tuple of:
          - ``sequence_tensor``  : np.ndarray of shape (N, 30, 4, 1)
          - ``aux_features``     : np.ndarray of shape (N, M)
          - ``guide_ids``        : List of N guide ID strings
    """
    df       = pd.read_csv(feature_csv)
    colnames = list(df.columns)

    # Locate the per-position nucleotide columns
    idx_A = [colnames.index(f"A{x}") for x in range(1, 31)]
    idx_T = [colnames.index(f"T{x}") for x in range(1, 31)]
    idx_G = [colnames.index(f"G{x}") for x in range(1, 31)]
    idx_C = [colnames.index(f"C{x}") for x in range(1, 31)]

    has_dinuc = "AA1" in colnames
    dinuc_pairs = [
        ("AA", "AT", "AG", "AC"),
        ("TA", "TT", "TG", "TC"),
        ("GA", "GT", "GG", "GC"),
        ("CA", "CT", "CG", "CC"),
    ]
    dinuc_indices: List[List[int]] = []
    if has_dinuc:
        for group in dinuc_pairs:
            for di in group:
                dinuc_indices.append(
                    [colnames.index(f"{di}{x}") for x in range(1, 30)]
                )

    id_index   = colnames.index("ID")
    all_seq_idx = (
        idx_A + idx_T + idx_G + idx_C
        + [i for group in dinuc_indices for i in group]
    )
    aux_indices = [
        i for i in range(len(colnames))
        if i not in all_seq_idx and i != id_index
    ]

    sequence_tensors: List[List] = []
    aux_rows:         List[List] = []
    guide_ids:        List[str]  = []

    for _, row in df.iterrows():
        guide_ids.append(str(row.iloc[id_index]))
        vals = row.tolist()

        a_vec = [int(vals[i]) for i in idx_A]
        t_vec = [int(vals[i]) for i in idx_T]
        g_vec = [int(vals[i]) for i in idx_G]
        c_vec = [int(vals[i]) for i in idx_C]

        onehot = [[a_vec[i], c_vec[i], g_vec[i], t_vec[i]] for i in range(30)]
        sequence_tensors.append(onehot)

        aux_row = [float(vals[i]) for i in aux_indices]
        if has_dinuc:
            kmer_flags = []
            seq_str = ""
            for i in range(30):
                if a_vec[i]:   seq_str += "A"
                elif t_vec[i]: seq_str += "T"
                elif c_vec[i]: seq_str += "C"
                elif g_vec[i]: seq_str += "G"
            for nuc1, nuc2 in itertools.product("ATCG", repeat=2):
                dinuc = nuc1 + nuc2
                for kmer in [seq_str[i:i+2] for i in range(len(seq_str)-1)]:
                    kmer_flags.append(1 if kmer == dinuc else 0)
            aux_row += kmer_flags

        aux_rows.append(aux_row)

    seq_array = np.asarray(sequence_tensors).reshape(-1, 30, 4, 1)
    aux_array = np.asarray(aux_rows)
    return seq_array, aux_array, guide_ids


# ---------------------------------------------------------------------------
# Step 11 — Run ensemble model and write scores
# ---------------------------------------------------------------------------

def run_model_ensemble(model_dir: str, seq_array: np.ndarray,
                       aux_array: np.ndarray) -> np.ndarray:
    """
    Load all ``.h5`` model files in *model_dir*, run inference, and return
    the ensemble mean prediction.

    Args:
        model_dir:  Directory containing trained Keras ``.h5`` model files.
        seq_array:  One-hot encoded sequence tensor (N, 30, 4, 1).
        aux_array:  Auxiliary feature matrix (N, M).

    Returns:
        1-D numpy array of length N with ensemble mean DGD scores.

    Raises:
        FileNotFoundError : If *model_dir* does not exist or contains no models.
    """
    model_files = [
        os.path.join(model_dir, f)
        for f in sorted(os.listdir(model_dir))
        if f.endswith(".h5")
    ]
    if not model_files:
        raise FileNotFoundError(
            f"No .h5 model files found in '{model_dir}'"
        )

    predictions: List[np.ndarray] = []
    for model_path in model_files:
        model = keras_models.load_model(model_path)
        pred  = model.predict([seq_array, aux_array])
        predictions.append(pred)

    ensemble_mean = np.mean(predictions, axis=0).flatten()
    logger.info(
        "Ensemble of %d models scored %d guides.", len(model_files), len(ensemble_mean)
    )
    return ensemble_mean


def score_guides(model_dir: str, output_path: str) -> None:
    """
    Load features, run the DGD ensemble model, merge scores with guide metadata,
    and write the final results CSV.

    Args:
        model_dir:   Path to the directory containing trained ``.h5`` model files.
        output_path: Destination path for the output CSV (default: ``DGD.csv``).

    Reads  : ``Deep_learning_file.csv``, ``Structure_file.csv``
    Writes : *output_path*
    """
    seq_array, aux_array, guide_ids = prepare_model_inputs(PipelineFiles.DEEP_LEARNING)
    scores = run_model_ensemble(model_dir, seq_array, aux_array)

    score_df   = pd.DataFrame({"ID": guide_ids, "DGD": scores})
    struct_df  = pd.read_csv(PipelineFiles.STRUCTURE_FILE)
    results    = struct_df.merge(score_df, on="ID", how="inner")
    results.to_csv(output_path, index=False)

    logger.info(
        "DGD scoring complete: %d guides written to '%s'",
        len(results), output_path,
    )


# ---------------------------------------------------------------------------
# Pipeline orchestrator
# ---------------------------------------------------------------------------

def run_pipeline(fasta_path: str, model_dir: str, output_path: str) -> None:
    """
    Run the full DGD pipeline from FASTA input to scored output CSV.

    Pipeline steps:
      1. Scan FASTA for guide candidates  → Structure_file.csv
      2. Build FASTA for RNAfold          → Structure_Connection.fa
      3. Compute target sequence features → Target_sequence_feature.csv
      4. RNAfold secondary structure      (external: RNAfold + b2ct)
      5. Parse connection matrix          → Structure_out.txt
      6. Run connection_to_matrix         (external C++ binary)
      7. Extract spacer–scaffold pairs    → spacer_scaffold_basepairs.csv
      8. Compute connection frequency     → spacer_scaffold_feature.csv
      9. Annotate structure regions       → Structural_annotation.csv
     10. Build connectivity features      → Feature_Data_Spacer_Scaffold.csv
     11. Compute final structural feats   → Deep_learning_file.csv
     12. Score with DGD ensemble          → DGD.csv

    Args:
        fasta_path:  Path to the input FASTA file.
        model_dir:   Path to the directory containing trained model ``.h5`` files.
        output_path: Destination path for final scores CSV.
    """
    logger.info("=== DGD Pipeline Starting ===")
    logger.info("Input:  %s", fasta_path)
    logger.info("Models: %s", model_dir)
    logger.info("Output: %s", output_path)

    # Step 1 — Scan FASTA for guide candidates
    logger.info("[1/12] Scanning guide candidates...")
    scan_guides(fasta_path)

    # Step 2 — Build FASTA for RNAfold
    logger.info("[2/12] Building RNAfold FASTA...")
    make_fasta_for_rnafold()

    # Step 3 — Compute target sequence features
    logger.info("[3/12] Computing target sequence features...")
    compute_target_features()

    # Step 4 — Run RNAfold (external binary)
    logger.info("[4/12] Running RNAfold...")
    _run_command(
        f"RNAfold -j0 --noPS < {PipelineFiles.FASTA_OUT} > {PipelineFiles.RNAFOLD_OUT}",
        shell=True,
        description="RNAfold",
    )
    _run_command(
        f"b2ct < {PipelineFiles.RNAFOLD_OUT} > {PipelineFiles.RNAFOLD_OUTS}",
        shell=True,
        description="b2ct",
    )

    # Step 5 — Parse connection matrix
    logger.info("[5/12] Parsing connection matrix...")
    parse_rnafold_output()

    # Step 6 — Run connection_to_matrix (compiled C++ binary)
    logger.info("[6/12] Running connection_to_matrix...")
    _run_command(
        f"./connection_to_matrix {PipelineFiles.STRUCTURE_OUT} 102 > {PipelineFiles.BASEPAIRS_CSV}",
        shell=True,
        description="connection_to_matrix",
    )

    # Step 7 — Extract spacer–scaffold pairs
    logger.info("[7/12] Extracting spacer–scaffold pairs...")
    extract_spacer_scaffold_pairs()

    # Step 8 — Compute connection frequency
    logger.info("[8/12] Computing connection frequency...")
    compute_connection_frequency()

    # Step 9 — Annotate structure regions
    logger.info("[9/12] Annotating structure regions...")
    annotate_structure_regions()

    # Step 10 — Build connectivity feature matrix
    logger.info("[10/12] Building connectivity features...")
    build_features()

    # Step 11 — Compute final structural features
    logger.info("[11/12] Computing final structural features...")
    compute_final_features()

    # Step 12 — Score with DGD ensemble
    logger.info("[12/12] Scoring with DGD model ensemble...")
    score_guides(model_dir, output_path)

    logger.info("=== DGD Pipeline Complete ===")
    logger.info("Results written to '%s'", output_path)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _run_command(cmd: str, shell: bool = False, description: str = "") -> None:
    """
    Run an external shell command via subprocess, raising on failure.

    Using ``subprocess.run`` (instead of ``os.system``) gives us:
      - Return code checking
      - Captured stderr for better error messages
      - No shell injection risk when ``shell=False``

    Args:
        cmd:         Command string or list of arguments.
        shell:       If True, run through the system shell (needed for pipes/redirects).
        description: Human-readable name for error messages.

    Raises:
        RuntimeError : If the command returns a non-zero exit code.
    """
    label = description or (cmd if isinstance(cmd, str) else " ".join(cmd))
    result = subprocess.run(
        cmd,
        shell=shell,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command '{label}' failed (exit {result.returncode}):\n{result.stderr}"
        )
    if result.stderr:
        logger.debug("stderr from '%s': %s", label, result.stderr.strip())


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser for DGD."""
    parser = argparse.ArgumentParser(
        prog="DGD",
        description=(
            "DGD — Deep Guide Designer\n"
            "CNN-based on-target scoring for CRISPR-Cas9 guide RNAs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python DGD.py input.fa\n"
            "  python DGD.py input.fa --output results.csv\n"
            "  python DGD.py input.fa --models ./my_models --verbose\n"
        ),
    )
    parser.add_argument(
        "fasta",
        metavar="FASTA",
        help="Input FASTA file with one or more sequences (100–10 000 nt each).",
    )
    parser.add_argument(
        "--output", "-o",
        metavar="CSV",
        default=PipelineFiles.OUTPUT,
        help=f"Output CSV file path (default: {PipelineFiles.OUTPUT}).",
    )
    parser.add_argument(
        "--models", "-m",
        metavar="DIR",
        default="./models",
        help="Directory containing trained DGD model .h5 files (default: ./models).",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose/debug logging.",
    )
    return parser


def main() -> None:
    """Entry point for the DGD command-line tool."""
    parser = _build_parser()
    args   = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        run_pipeline(
            fasta_path=args.fasta,
            model_dir=args.models,
            output_path=args.output,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as err:
        logger.error("%s", err)
        sys.exit(1)


if __name__ == "__main__":
    main()
