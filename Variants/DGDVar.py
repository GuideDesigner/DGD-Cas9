#!/usr/bin/env python3
"""
Variants/DGDVar.py — DGD pipeline for broad-PAM Cas9 variant scoring
=====================================================================
Authors : Vipin Menon, Jang-il Sohn, Seokju Park, Jin-Wu Nam
Lab     : Bioinformatics & Genomics Lab, Hanyang University, Seoul 04763, Korea
Contact : a.vipin.menon@gmail.com | jwnam@hanyang.ac.kr
Original: August 2021 | Modernized: 2026

Scores CRISPR-Cas9 sgRNA candidates against 9 SpCas9 variant models using the
DGD CNN ensemble. Unlike DGD.py (NGG-only), DGDVar scans ALL 16 dinucleotide
PAMs to support xCas9, SpCas9-NG, and other broad-PAM Cas9 variants.

Usage:
    python DGDVar.py input.fa
    python DGDVar.py input.fa --output results.csv --models ./models

Output columns (DGDVar.csv):
    DGDSpCas9, DGDeSpCas9, DGDHypaCas9, DGDSpCas9Hf1, DGDSniperCas9,
    DGDevoCas9, DGDxCas9, DGDSpCas9VRQR, DGDSpCas9NG
"""

import argparse
import logging
import os
import subprocess
import sys
from collections import OrderedDict, defaultdict
from itertools import chain
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import RNA
from Bio.SeqUtils import MeltingTemp as mt
from tensorflow.keras import models as tf_models

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from make_arrays import (
    adjust_partner_index,
    return_connection_length,
    return_dimers_array,
    return_masked_array,
    return_partner_base,
    return_partner_index,
)
from sequence_utils import parse_fasta, reverse_complement
from stacking_model import return_stacking_model

logger = logging.getLogger(__name__)

# Scaffold sequence for SpCas9 variants (xCas9/SpCas9-NG optimized)
_SCAFFOLD_VAR: str = (
    "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCC"
    "GTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
)

# All 16 dinucleotide PAM sequences for broad-PAM scanning
_ALL_PAMS: List[str] = [
    "GG", "AA", "AT", "AG", "AC", "GA", "GT", "GC",
    "TA", "TT", "TG", "TC", "CA", "CT", "CG", "CC",
]


class PipelineFiles:
    """Centralised store for all intermediate pipeline file paths."""
    STRUCTURE_FILE:             str = "Structure_file.csv"
    STRUCTURE_CONN_FA:          str = "Structure_Connection.fa"
    STRUCTURE_CONN_CSV:         str = "Structure_Connection.csv"
    STRUCTURE_CONN_OUT:         str = "Structure_Connection.out"
    STRUCTURE_CONN_OUTS:        str = "Structure_Connection.outs"
    STRUCTURE_OUT:              str = "Structure_out.txt"
    STRUCTURE_BASEPAIRS:        str = "Structure_basepairs.csv"
    SPACER_SCAFFOLD_BP:         str = "spacer_scaffold_basepairs.csv"
    SPACER_SCAFFOLD_FEATURE:    str = "spacer_scaffold_feature.csv"
    STRUCTURAL_ANNOTATION:      str = "Structural_annotation.csv"
    TARGET_FEATURES:            str = "Target_sequence_feature.csv"
    FEATURE_DATA:               str = "Feature_Data_Spacer_Scaffold.csv"
    DEEP_LEARNING_FILE:         str = "Deep_learning_file.csv"
    STRUCTURE_CAS9_TMP:         str = "Structure_Cas9_out.txt"


def _run_command(cmd: List[str], *, shell: bool = False, description: str = "") -> None:
    """Run an external command, raising RuntimeError on failure."""
    result = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
    if result.returncode != 0:
        msg = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            f"Command failed{' (' + description + ')' if description else ''}: {msg}"
        )


def scan_guides(fasta_path: str, output_path: str = PipelineFiles.STRUCTURE_FILE) -> None:
    """
    Step 1: Scan a FASTA file for all guide candidates across all 16 dinucleotide PAMs.

    Args:
        fasta_path: Input FASTA file (100–10,000 nt sequences).
        output_path: Output CSV path for guide candidates.
    """
    guide_fwd: Dict[str, list] = {}
    guide_rev: Dict[str, list] = {}

    for seq_id, sequence in parse_fasta(fasta_path):
        seq_len = len(sequence)
        if not (100 <= seq_len <= 10_000):
            logger.warning("Skipping %s: length %d outside [100, 10000]", seq_id, seq_len)
            continue

        for pam in _ALL_PAMS:
            rc_pam = reverse_complement(pam)
            pos = 0
            while True:
                st = sequence.find(pam, pos)
                if st == -1:
                    break
                pos = st + 1
                if 25 < st < seq_len - 5:
                    nstart = st - 25
                    nend   = st + 5
                    seq_win = sequence[nstart:nend]
                    idn = f"{seq_id}:{nstart}:{nend}:+"
                    guide_fwd[idn] = [nstart, nend, "+", seq_win]

            pos = 0
            while True:
                st = sequence.find(rc_pam, pos)
                if st == -1:
                    break
                pos = st + 1
                if 3 < st < seq_len - 27:
                    nstart = st - 3
                    nend   = st + 27
                    seq_win = sequence[nstart:nend]
                    rc_win  = reverse_complement(seq_win)
                    idn = f"{seq_id}:{nstart}:{nend}:-"
                    guide_rev[idn] = [nstart, nend, "-", rc_win]

    all_guides = dict(chain(guide_fwd.items(), guide_rev.items()))
    sorted_guides = OrderedDict(sorted(all_guides.items(), key=lambda x: x[1][3], reverse=True))

    with open(output_path, "w", encoding="utf-8") as out:
        out.write("ID,Start,End,Strand,Sequence\n")
        for gid, vals in sorted_guides.items():
            out.write(f"{gid},{vals[0]},{vals[1]},{vals[2]},{vals[3]}\n")

    logger.info("Wrote %d guide candidates → %s", len(sorted_guides), output_path)


def make_fasta_for_rnafold(
    structure_file: str = PipelineFiles.STRUCTURE_FILE,
    fasta_out: str = PipelineFiles.STRUCTURE_CONN_FA,
    csv_out: str = PipelineFiles.STRUCTURE_CONN_CSV,
) -> None:
    """
    Step 2: Append the sgRNA scaffold to each guide (positions 4:24) and write FASTA.

    Args:
        structure_file: Guide candidate CSV.
        fasta_out: Output FASTA file for RNAfold.
        csv_out: Output CSV with ID and full sequence.
    """
    with (
        open(structure_file, "r", encoding="utf-8") as fh,
        open(fasta_out, "w", encoding="utf-8") as fa,
        open(csv_out, "w", encoding="utf-8") as csv_fh,
    ):
        csv_fh.write("ID,Sequence\n")
        lines = fh.readlines()[1:]  # skip header
        for line in lines:
            info = line.strip().split(",")
            seq_id = info[0]
            guide  = info[4][4:24] + _SCAFFOLD_VAR
            fa.write(f">{seq_id}\n{guide}\n")
            csv_fh.write(f"{seq_id},{guide}\n")

    logger.info("Wrote scaffold-appended sequences → %s", fasta_out)


def compute_target_features(
    structure_file: str = PipelineFiles.STRUCTURE_FILE,
    output_path: str = PipelineFiles.TARGET_FEATURES,
) -> None:
    """
    Step 3: Compute one-hot encoding, entropy, free energy, GC content, and Tm for each guide.

    Args:
        structure_file: Guide candidate CSV.
        output_path: Output feature CSV.
    """
    nuc = ["A", "T", "G", "C"]
    dinuct: List[str] = [f"{b}{a}" for a in nuc for b in nuc]

    header = (
        [f"A{i}" for i in range(1, 31)] +
        [f"T{i}" for i in range(1, 31)] +
        [f"G{i}" for i in range(1, 31)] +
        [f"C{i}" for i in range(1, 31)] +
        [f"{d}{i}" for d in ["AA","TA","GA","CA","AT","TT","GT","CT",
                              "AG","TG","GG","CG","AC","TC","GC","CC"]
         for i in range(1, 30)]
    )
    extra_cols = (
        "Entropy,Energy,GCcount,Gchigh,GClow,MeltingTemperature,"
        "A,T,G,C,AA,AT,AG,AC,CA,CG,CC,CT,GA,GC,GG,GT,TA,TC,TG,TT"
    )

    merged: defaultdict = defaultdict(list)
    sigma: defaultdict = defaultdict(list)

    nt: Dict[str, list] = {}
    dt: Dict[str, list] = {}
    ent: Dict[str, str] = {}
    ene: Dict[str, str] = {}
    gc: Dict[str, str] = {}
    gc_high: Dict[str, str] = {}
    gc_low: Dict[str, str] = {}
    Tm_map: Dict[str, str] = {}
    PosA: Dict[str, str] = {}
    PosT: Dict[str, str] = {}
    PosG: Dict[str, str] = {}
    PosC: Dict[str, str] = {}
    PosAA: Dict[str, str] = {}
    PosAT: Dict[str, str] = {}
    PosAG: Dict[str, str] = {}
    PosAC: Dict[str, str] = {}
    PosCA: Dict[str, str] = {}
    PosCC: Dict[str, str] = {}
    PosCG: Dict[str, str] = {}
    PosCT: Dict[str, str] = {}
    PosGA: Dict[str, str] = {}
    PosGC: Dict[str, str] = {}
    PosGG: Dict[str, str] = {}
    PosGT: Dict[str, str] = {}
    PosTA: Dict[str, str] = {}
    PosTC: Dict[str, str] = {}
    PosTG: Dict[str, str] = {}
    PosTT: Dict[str, str] = {}

    with open(structure_file, "r", encoding="utf-8") as fh:
        lines = fh.readlines()[1:]

    for line in lines:
        seq = line.strip().split(",")
        ids              = seq[0]
        complete_seq     = seq[4]
        target_seq       = complete_seq[4:24].upper()
        extended_seq     = complete_seq.upper()

        nt[ids] = [str(int(nuc[y] == extended_seq[x]))
                   for y in range(4) for x in range(len(extended_seq))]
        dt[ids] = [str(int(dinuct[w] == extended_seq[z:z+2]))
                   for w in range(16) for z in range(len(extended_seq)-1)]

        nts = {n: extended_seq.count(n) / len(target_seq) for n in nuc}
        ent[ids] = str(round(sum(
            -(p * np.log2(p)) for p in nts.values() if p != 0
        ), 1))
        ene[ids] = str(round(RNA.fold(target_seq)[-1], 0))

        PosA[ids]  = str(extended_seq.count("A"))
        PosT[ids]  = str(extended_seq.count("T"))
        PosG[ids]  = str(extended_seq.count("G"))
        PosC[ids]  = str(extended_seq.count("C"))
        Tm_map[ids] = str(mt.Tm_NN(target_seq))

        for dinuc, dmap in [
            ("AA", PosAA), ("AT", PosAT), ("AG", PosAG), ("AC", PosAC),
            ("CA", PosCA), ("CC", PosCC), ("CG", PosCG), ("CT", PosCT),
            ("GA", PosGA), ("GC", PosGC), ("GG", PosGG), ("GT", PosGT),
            ("TA", PosTA), ("TC", PosTC), ("TG", PosTG), ("TT", PosTT),
        ]:
            dmap[ids] = str(extended_seq.count(dinuc))

        gc_count = target_seq.count("G") + target_seq.count("C")
        gc[ids]      = str(round(gc_count / len(target_seq) * 100, 0))
        gc_low[ids]  = str(int(gc_count < 10))
        gc_high[ids] = str(int(gc_count >= 10))

    for d in [nt, dt]:
        for k in d:
            d[k] = ",".join(d[k])

    dict_list = [
        nt, dt, ent, ene, gc, gc_high, gc_low, Tm_map,
        PosA, PosT, PosG, PosC,
        PosAA, PosAT, PosAG, PosAC,
        PosCA, PosCC, PosCG, PosCT,
        PosGA, PosGC, PosGG, PosGT,
        PosTA, PosTC, PosTG, PosTT,
    ]
    for d in dict_list:
        for k, v in d.items():
            merged[k].append(v)
    for k, v in merged.items():
        sigma[k] = ",".join(v)

    with open(output_path, "w", encoding="utf-8") as out:
        out.write("ID," + ",".join(header) + "," + extra_cols + "\n")
        for k, v in sorted(sigma.items()):
            out.write(f"{k},{v}\n")

    logger.info("Wrote target features → %s", output_path)


def parse_rnafold_output(
    input_path: str = PipelineFiles.STRUCTURE_CONN_OUTS,
    output_path: str = PipelineFiles.STRUCTURE_OUT,
) -> None:
    """
    Step 5: Parse b2ct-formatted RNAfold output into a tabular connection matrix.

    Args:
        input_path: b2ct output file (Structure_Connection.outs).
        output_path: Output connection matrix (Structure_out.txt).
    """
    header = ["ID"] + [f"Pos{i}" for i in range(1, 113)]
    sara: Dict[str, list] = {}
    current_key: str = ""

    with open(input_path, "r", encoding="utf-8") as fh:
        for line in fh:
            info = list(filter(None, line.strip().split()))
            if len(info) == 5:
                current_key = info[4]
                sara[current_key] = []
            elif current_key:
                sara[current_key].append(info[4])

    tmp_path = PipelineFiles.STRUCTURE_CAS9_TMP
    df = pd.DataFrame.from_dict(sara, orient="index")
    df.to_csv(tmp_path, sep="\t", header=False)
    df = pd.read_csv(tmp_path, sep="\t", names=header)
    df.to_csv(output_path, sep="\t", index=False)
    os.remove(tmp_path)
    logger.info("Wrote connection matrix → %s", output_path)


def extract_spacer_scaffold_pairs(
    input_path: str = PipelineFiles.STRUCTURE_BASEPAIRS,
    output_path: str = PipelineFiles.SPACER_SCAFFOLD_BP,
) -> None:
    """
    Step 7: Extract spacer–scaffold base-pair columns from the full matrix.

    Args:
        input_path: Full connection matrix CSV.
        output_path: Output spacer–scaffold base-pair CSV.
    """
    conn_d = pd.read_csv(input_path)
    columns = [f"Connection_Pos{i}_Pos{j}"
               for i in range(1, 21) for j in range(21, 113)] + ["ID"]
    pd.DataFrame(conn_d, columns=columns).to_csv(output_path, index=False)
    logger.info("Wrote spacer–scaffold pairs → %s", output_path)


def compute_connection_frequency(
    input_path: str = PipelineFiles.SPACER_SCAFFOLD_BP,
    output_path: str = PipelineFiles.SPACER_SCAFFOLD_FEATURE,
) -> None:
    """
    Step 8: Build per-position connection frequency metadata.

    Args:
        input_path: Spacer–scaffold base-pair CSV.
        output_path: Output position-frequency CSV.
    """
    connection = pd.read_csv(input_path)
    cons = connection.set_index("ID").T.reset_index()
    cons_info = cons[[cons.columns[0]]].copy()
    cons_info.columns = ["nucleotide"]
    split = pd.DataFrame(
        [x.split("_") for x in cons_info["nucleotide"].tolist()]
    )
    cons_info["Pos_A"] = split.iloc[:, 1].str.extract(r"(\d+)", expand=False)
    cons_info["Pos_B"] = split.iloc[:, 2].str.extract(r"(\d+)", expand=False)
    cons_info.to_csv(output_path, index=False)
    logger.info("Wrote connection frequency → %s", output_path)


def annotate_structure_regions(
    input_path: str = PipelineFiles.SPACER_SCAFFOLD_FEATURE,
    output_path: str = PipelineFiles.STRUCTURAL_ANNOTATION,
) -> None:
    """
    Step 9: Annotate each position with its structural region (R, TL, AR, LR, SL1, SL2, SL3, NS).

    Args:
        input_path: Connection frequency CSV.
        output_path: Annotated output CSV.
    """
    df = pd.read_csv(input_path)
    df["Pos_B"] = pd.to_numeric(df["Pos_B"])

    region_ranges = [
        ("TL",  33,  36),
        ("SL1", 54,  58),
        ("SL2", 73,  76),
        ("SL3", 88,  90),
        ("R",   21,  32),
        ("AR",  37,  49),
        ("LR",  63,  67),
    ]
    for label, lo, hi in region_ranges:
        df.loc[df["Pos_B"].between(lo, hi), "Structure"] = label
    df["Structure"].fillna("NS", inplace=True)
    df.to_csv(output_path, index=False)
    logger.info("Wrote structural annotation → %s", output_path)


def _build_structure_df(
    region_df: pd.DataFrame,
    connection_bp: pd.DataFrame,
    prefix: str,
) -> pd.DataFrame:
    """
    Build a binary connectivity DataFrame for one structural region.

    For each unique spacer position in `region_df`, determine which guide IDs
    have at least one base-pair to the scaffold in that region.

    Args:
        region_df:     Subset of annotation DataFrame for one structural region.
        connection_bp: Full spacer–scaffold base-pair matrix.
        prefix:        Column-name prefix (e.g. "AR", "TL").

    Returns:
        DataFrame with one column per spacer position plus a ``Total_{prefix}`` column.
    """
    positions = sorted(set(region_df["Pos_A"].tolist()), reverse=True)
    nested: Dict[str, Dict[str, int]] = {}

    for pos in positions:
        col_name  = f"{prefix}{pos}"
        nuc_list  = list(set(
            region_df[region_df["Pos_A"] == pos]["nucleotide"].tolist()
        ))
        conn_info = connection_bp[nuc_list].copy()
        conn_info.insert(0, "ID", connection_bp[["ID"]])
        id_dict   = conn_info.set_index("ID").T.to_dict("list")
        nested[col_name] = {k: int(sum(v) > 0) for k, v in id_dict.items()}

    result_df = pd.DataFrame(nested)
    result_df[f"Total_{prefix}"] = result_df.sum(axis=1)
    result_df.index.name = "ID"
    result_df.reset_index(inplace=True)
    return result_df


def build_features(
    bp_input: str = PipelineFiles.SPACER_SCAFFOLD_BP,
    seq_input: str = PipelineFiles.TARGET_FEATURES,
    annot_input: str = PipelineFiles.STRUCTURAL_ANNOTATION,
    output_path: str = PipelineFiles.FEATURE_DATA,
) -> None:
    """
    Step 10: Build per-region structural connectivity features and merge with sequence features.

    Args:
        bp_input:    Spacer–scaffold base-pair CSV.
        seq_input:   Target sequence features CSV.
        annot_input: Structural annotation CSV.
        output_path: Output combined feature CSV.
    """
    conn_bp  = pd.read_csv(bp_input)
    seq_df   = pd.read_csv(seq_input)
    annot_df = pd.read_csv(annot_input)

    regions = ["AR", "NS", "R", "SL1", "LR", "SL2", "SL3", "TL"]
    region_dfs: Dict[str, pd.DataFrame] = {}

    for region in regions:
        sub = annot_df[annot_df["Structure"] == region][
            ["nucleotide", "Pos_A", "Pos_B", "Structure"]
        ]
        region_dfs[region] = _build_structure_df(sub, conn_bp, region)

    conn_totals = pd.concat(
        [region_dfs[r][f"Total_{r}"] for r in regions], axis=1, join="outer"
    )
    total = conn_totals.sum(axis=1).rename("Total_Connection")

    ar = region_dfs["AR"]
    drop_cols = ["ID"]
    merged_conn = pd.concat(
        [ar] + [region_dfs[r].drop(columns=drop_cols) for r in regions if r != "AR"] + [total],
        axis=1, join="outer"
    )
    seq_df.merge(merged_conn, on="ID", how="inner").to_csv(output_path, index=False)
    logger.info("Wrote combined features → %s", output_path)


def _read_dot_bracket_file(file_path: str) -> Dict[str, list]:
    """
    Parse an RNAfold dot-bracket format file into a dict.

    Returns:
        {sequence_id: [sequence, structure_without_energy, energy_string]}
    """
    result: Dict[str, list] = {}
    with open(file_path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()

    for i, line in enumerate(lines):
        if i % 3 == 0:
            seq_id = line[1:].strip()
            result[seq_id] = []
        elif i % 3 == 1:
            result[seq_id].append(line.rstrip("\n"))
        elif i % 3 == 2:
            result[seq_id].append(line[:line.rfind("(")])
            result[seq_id].append(line[line.rfind("("):line.find("\n")])
    return result


def _count_monomer(sequence: str, partner_index: list) -> Dict[str, int]:
    monomer = {"A": 0, "C": 0, "G": 0, "U": 0}
    for i in range(20):
        if partner_index[i] != 0:
            monomer[sequence[i]] += 1
            monomer[sequence[partner_index[i]]] += 1
    return monomer


def _count_dimer(dimers: list) -> Dict[str, int]:
    dinucleotide = {
        "AA": 0, "AC": 0, "AG": 0, "AU": 0,
        "CA": 0, "CC": 0, "CG": 0, "CU": 0,
        "GA": 0, "GC": 0, "GG": 0, "GU": 0,
        "UA": 0, "UC": 0, "UG": 0, "UU": 0,
    }
    for d in dimers:
        dinucleotide[d] += 1
    return dinucleotide


def _count_gc_content(sequence: str, partner_index: list) -> int:
    count = 0
    for i in range(20):
        if partner_index[i] != 0:
            pair = {sequence[i], sequence[partner_index[i]]}
            if pair == {"G", "C"}:
                count += 1
    return count


def _calculate_free_energy(dimers: list) -> float:
    stacking = return_stacking_model()
    energy = 0.0
    for i in range(0, len(dimers), 2):
        energy += stacking[dimers[i]][dimers[i + 1][::-1]]
    return energy


def compute_final_features(
    rnafold_out: str = PipelineFiles.STRUCTURE_CONN_OUT,
    feature_input: str = PipelineFiles.FEATURE_DATA,
    output_path: str = PipelineFiles.DEEP_LEARNING_FILE,
) -> None:
    """
    Step 11: Compute spacer–scaffold structural features from dot-bracket and merge.

    Args:
        rnafold_out:   RNAfold dot-bracket output file.
        feature_input: Combined feature CSV.
        output_path:   Output deep-learning input CSV.
    """
    id_data    = _read_dot_bracket_file(rnafold_out)
    feat_df    = pd.read_csv(feature_input)

    ids_gc: Dict[str, list] = {}
    ids_mono_di: Dict[str, list] = {}

    for seq_id, vals in id_data.items():
        sequence, structure = vals[0], vals[1]
        p_idx    = adjust_partner_index(return_partner_index(structure))
        p_base   = return_partner_base(p_idx, sequence)
        _, count = return_connection_length(p_idx)
        dimers   = return_dimers_array(sequence, p_idx)

        gc  = _count_gc_content(sequence, p_idx) / count if count else 0.0
        eng = _calculate_free_energy(dimers)
        ids_gc[seq_id]      = [gc, eng]
        ids_mono_di[seq_id] = [_count_monomer(sequence, p_idx), _count_dimer(dimers)]

    mono_di = {
        k: {key: val for d in v for key, val in d.items()}
        for k, v in ids_mono_di.items()
    }
    df1 = pd.DataFrame(mono_di).T
    df1.index.name = "ID"
    df1.reset_index(inplace=True)

    df2 = pd.DataFrame.from_dict(ids_gc).T
    df2.columns = ["GC_ratio", "Gibbs_Energy"]
    df2.index.name = "ID"
    df2.reset_index(inplace=True)

    combined = pd.merge(df1, df2, on="ID", how="inner")
    combined.columns = [
        "ID",
        "Spacer_Scaffold_A", "Spacer_Scaffold_C", "Spacer_Scaffold_G", "Spacer_Scaffold_U",
        "Spacer_Scaffold_AA", "Spacer_Scaffold_AC", "Spacer_Scaffold_AG", "Spacer_Scaffold_AU",
        "Spacer_Scaffold_CA", "Spacer_Scaffold_CC", "Spacer_Scaffold_CG", "Spacer_Scaffold_CU",
        "Spacer_Scaffold_GA", "Spacer_Scaffold_GC", "Spacer_Scaffold_GG", "Spacer_Scaffold_GU",
        "Spacer_Scaffold_UA", "Spacer_Scaffold_UC", "Spacer_Scaffold_UG", "Spacer_Scaffold_UU",
        "Spacer_Scaffold_GC_ratio", "Spacer_Scaffold_Gibbs_Energy",
    ]
    feat_df.merge(combined, on="ID").to_csv(output_path, index=False)
    logger.info("Wrote final features → %s", output_path)


def _get_input(input_path: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Parse Deep_learning_file.csv into model input arrays.

    Returns:
        Tuple of (one-hot sequence array, additional feature array, guide IDs)
    """
    import itertools

    with open(input_path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()

    colnames    = [x.strip('"') for x in lines[0].strip().split(",")]
    data_lines  = lines[1:]
    idss: List[str] = [l.strip().split(",")[0] for l in data_lines]

    idx_a = [colnames.index(f"A{x}") for x in range(1, 31)]
    idx_t = [colnames.index(f"T{x}") for x in range(1, 31)]
    idx_g = [colnames.index(f"G{x}") for x in range(1, 31)]
    idx_c = [colnames.index(f"C{x}") for x in range(1, 31)]

    idx_di: List[int] = []
    if "AA1" in colnames:
        for pair in ["AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA","CT","CG","CC"]:
            idx_di += [colnames.index(f"{pair}{x}") for x in range(1, 30)]

    idx_id     = [colnames.index("ID")]
    idx_others = [i for i in range(len(colnames))
                  if i not in idx_a + idx_t + idx_g + idx_c + idx_id + idx_di]

    in_: List[list] = []
    add_: List[list] = []

    for line in data_lines:
        row = line.strip().split(",")
        _a = [int(row[i]) for i in idx_a]
        _t = [int(row[i]) for i in idx_t]
        _g = [int(row[i]) for i in idx_g]
        _c = [int(row[i]) for i in idx_c]
        _oth = [float(row[i]) for i in idx_others]

        seq = "".join(
            "A" if _a[i] else "T" if _t[i] else "C" if _c[i] else "G" if _g[i] else "?"
            for i in range(len(_a))
        )
        onehot = [[_a[i], _c[i], _g[i], _t[i]] for i in range(len(_a))]
        dinuc = [
            int("".join(nuc) == seq[j:j+2])
            for nuc in itertools.product("ATCG", repeat=2)
            for j in range(len(seq)-1)
        ]
        _oth += dinuc
        in_.append(onehot)
        add_.append(_oth)

    in_arr  = np.asarray(in_).reshape(-1, 30, 4, 1)
    add_arr = np.asarray(add_)
    return in_arr, add_arr, idss


def _test_model(model_dir: str, x: np.ndarray, a: np.ndarray) -> np.ndarray:
    """Load all .h5 models from a directory and return their ensemble mean prediction."""
    results = []
    for fname in filter(lambda f: f.endswith(".h5"), os.listdir(model_dir)):
        model = tf_models.load_model(os.path.join(model_dir, fname))
        results.append(model.predict([x, a]))
    return np.mean(results, axis=0)


def score_guides(
    models_base_dir: str,
    input_path: str = PipelineFiles.DEEP_LEARNING_FILE,
    struct_input: str = PipelineFiles.STRUCTURE_FILE,
    output_path: str = "DGDVar.csv",
) -> None:
    """
    Step 12: Score guide candidates using 9 SpCas9 variant CNN ensembles.

    Args:
        models_base_dir: Base directory containing variant model subdirs.
        input_path:      Deep learning input CSV.
        struct_input:    Guide structure CSV.
        output_path:     Output scores CSV.
    """
    struct_df = pd.read_csv(struct_input)
    x, a, pop_id = _get_input(input_path)

    variant_dirs = {
        "DGDSpCas9":    os.path.join(models_base_dir, "SpCas9",    "models"),
        "DGDeSpCas9":   os.path.join(models_base_dir, "eSpCas9",   "models"),
        "DGDHypaCas9":  os.path.join(models_base_dir, "HypaCas9",  "models"),
        "DGDSpCas9Hf1": os.path.join(models_base_dir, "SpCas9-Hf1","models"),
        "DGDSniperCas9":os.path.join(models_base_dir, "Sniper-Cas9","models"),
        "DGDevoCas9":   os.path.join(models_base_dir, "evoCas9",   "models"),
        "DGDxCas9":     os.path.join(models_base_dir, "xCas9",     "models"),
        "DGDSpCas9VRQR":os.path.join(models_base_dir, "SpCas9-VRQR","models"),
        "DGDSpCas9NG":  os.path.join(models_base_dir, "SpCas9-NG", "models"),
    }

    score_frames: List[pd.DataFrame] = [pd.DataFrame(pop_id, columns=["ID"])]
    for col_name, model_dir in variant_dirs.items():
        result = _test_model(model_dir, x, a)
        score_frames.append(
            pd.DataFrame(result.reshape(len(result), 1), columns=[col_name])
        )

    score_df = pd.concat(score_frames, axis=1, join="outer")
    merged   = struct_df.merge(score_df, on="ID", how="inner")
    merged.to_csv(output_path, index=False)

    # Re-read and apply no_activity threshold
    dgd = pd.read_csv(output_path)
    for col in variant_dirs:
        dgd.loc[dgd[col] < 0.0, col] = "no_activity"
    dgd.to_csv(output_path, index=False)
    logger.info("Wrote variant scores → %s", output_path)


def run_pipeline(
    fasta_path: str,
    output_path: str = "DGDVar.csv",
    models_dir: str = ".",
) -> None:
    """
    Run the complete 12-step DGDVar pipeline.

    Args:
        fasta_path:  Input FASTA file.
        output_path: Output CSV path for variant scores.
        models_dir:  Directory containing variant model subdirectories.
    """
    steps = [
        ("Step  1 — Scanning guides",          lambda: scan_guides(fasta_path)),
        ("Step  2 — Building scaffold FASTA",   make_fasta_for_rnafold),
        ("Step  3 — Computing target features", compute_target_features),
        ("Step  4 — RNAfold",                   lambda: _run_command(
            "RNAfold -j0 --noPS <Structure_Connection.fa >Structure_Connection.out",
            shell=True, description="RNAfold"
        )),
        ("Step  5 — Parsing RNAfold output",    parse_rnafold_output),
        ("Step  6 — connection_to_matrix",      lambda: _run_command(
            "./connection_to_matrix Structure_out.txt 112 > Structure_basepairs.csv",
            shell=True, description="connection_to_matrix"
        )),
        ("Step  7 — Extracting spacer–scaffold pairs",    extract_spacer_scaffold_pairs),
        ("Step  8 — Computing connection frequency",      compute_connection_frequency),
        ("Step  9 — Annotating structural regions",       annotate_structure_regions),
        ("Step 10 — Building features",                   build_features),
        ("Step 11 — Computing final features",            compute_final_features),
        ("Step 12 — Scoring with variant ensembles",      lambda: score_guides(
            models_dir, output_path=output_path
        )),
    ]

    for label, step_fn in steps:
        logger.info(label)
        step_fn()

    logger.info("Pipeline complete. Results: %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "DGDVar — DGD pipeline for broad-PAM Cas9 variant scoring. "
            "Scores guides against 9 SpCas9 variant models."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("FASTA", help="Input FASTA file (sequences must be 100–10,000 nt)")
    parser.add_argument("--output", "-o", default="DGDVar.csv", help="Output CSV file")
    parser.add_argument("--models", "-m", default=".", help="Base directory for variant model subdirs")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )

    try:
        run_pipeline(args.FASTA, output_path=args.output, models_dir=args.models)
    except (FileNotFoundError, ValueError, RuntimeError) as err:
        logger.error("%s", err)
        sys.exit(1)


if __name__ == "__main__":
    main()
