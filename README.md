# DGD-Cas9 — Deep Guide Designer for CRISPR-Cas9

[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![TensorFlow 2.x](https://img.shields.io/badge/TensorFlow-2.x-orange.svg)](https://www.tensorflow.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**DGD-Cas9** is a deep-learning pipeline for designing and scoring CRISPR-Cas9 sgRNA guide candidates. It integrates sequence features, RNA secondary structure, and spacer–scaffold base-pairing connectivity into a CNN ensemble model to predict guide activity.

**Authors:** Vipin Menon, Jang-il Sohn, Seokju Park, Jin-Wu Nam
**Lab:** Bioinformatics & Genomics Lab (BIG Lab), Hanyang University, Seoul 04763, Korea
**Contact:** a.vipin.menon@gmail.com | jwnam@hanyang.ac.kr
**Original:** August 2021 | **Modernized:** 2026

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Variants](#pipeline-variants)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Details (12 Steps)](#pipeline-details-12-steps)
- [Variants — DGDVar.py](#variants--dgdvarpy-broad-pam-cas9-variants)
- [Base Editors — DGDbaseeditor.py](#base-editors--dgdbaseeditorpy)
- [Accessory Scripts](#accessory-scripts)
- [Output Files](#output-files)
- [Git Push Commands](#git-push-commands)
- [Citation](#citation)

---

## Overview

DGD-Cas9 predicts CRISPR-Cas9 sgRNA on-target activity using a CNN ensemble that integrates three complementary feature sets:

| Feature Group | Description |
|---|---|
| Target sequence | One-hot encoding, Shannon entropy, RNA free energy, GC content, melting temperature |
| Spacer–scaffold contacts | Binary connectivity across 8 structural scaffold regions (R, TL, AR, LR, SL1, SL2, SL3, NS) |
| Sequence context | Dinucleotide composition, Gibbs free energy, GC ratio of spacer–scaffold interface |

---

## Pipeline Variants

Three pipeline scripts are provided for different CRISPR applications:

| Script | Application | PAM | Guide window | Models | Output file |
|---|---|---|---|---|---|
| `DGD.py` | Standard SpCas9 | NGG | 30 nt | SpCas9 DGD ensemble | `DGD.csv` |
| `Variants/DGDVar.py` | Broad-PAM Cas9 variants | All 16 dinucleotides | 30 nt | 9 SpCas9 variant ensembles | `DGDVar.csv` |
| `Base-Editors (BE)/DGDbaseeditor.py` | Base editing (ABE/CBE) | NGG | 24 nt | ABE + CBE ensembles | `DGDBE.csv` |

---

## Repository Structure

```
DGD-Cas9/
│
├── DGD.py                          # Main pipeline — SpCas9 (NGG, 30 nt window)
├── sequence_utils.py               # Shared utilities: parse_fasta(), reverse_complement()
├── make_arrays.py                  # RNA structure array helpers
├── stacking_model.py               # RNA dinucleotide stacking free-energy model
├── connection_to_matrix            # Compiled C++ binary (run `make` to build)
├── connection_to_matrix.cpp        # C++ source for base-pair matrix builder
├── Makefile                        # Build file for C++ binary
├── requirements.txt
├── README.md
├── .gitignore
│
├── models/                         # Trained SpCas9 DGD .h5 model files
│
├── Accessory/                      # Standalone wrappers for each pipeline step
│   ├── get_sequence.py             # Reverse complement / FASTA utilities
│   ├── fastamaker.py               # Step  2: scaffold-appended FASTA
│   ├── targetsequence.py           # Step  3: target sequence features
│   ├── Connectstructure.py         # Step  5: parse RNAfold b2ct output
│   ├── spacerscaffold.py           # Step  7: spacer–scaffold base pairs
│   ├── spacerconnectionfrequency.py# Step  8: connection frequency table
│   ├── serialconnection.py         # Step  9: annotate structural regions
│   ├── featuremaker.py             # Step 10: structural feature matrix
│   ├── finalfeatures.py            # Step 11: final spacer–scaffold features
│   └── score_deep.py               # Step 12: DGD ensemble scoring
│
├── Variants/
│   └── DGDVar.py                   # Broad-PAM pipeline — 9 SpCas9 variant models
│
└── Base-Editors (BE)/
    └── DGDbaseeditor.py            # Base editor pipeline — ABE + CBE (24 nt window)
```

---

## Installation

### Prerequisites

- Python 3.8 or higher
- GCC (for compiling the C++ binary)
- RNAfold and b2ct from the ViennaRNA package

### 1. Clone the repository

```bash
git clone https://github.com/GuideDesigner/DGD-Cas9.git
cd DGD-Cas9
```

### 2. Create a virtual environment (recommended)

```bash
python3 -m venv venv
source venv/bin/activate        # macOS / Linux
venv\Scripts\activate           # Windows
```

### 3. Install Python dependencies

```bash
pip install -r requirements.txt
```

**`requirements.txt`:**
```
numpy>=1.21
pandas>=1.3
biopython>=1.79
tensorflow>=2.6
ViennaRNA>=2.5.0
```

ViennaRNA is best installed via conda if the pip build fails on your system:

```bash
conda install -c bioconda viennarna
```

### 4. Compile the C++ binary

```bash
make
```

This compiles `connection_to_matrix.cpp` → the `connection_to_matrix` executable used in Step 6 of the pipeline.

### 5. Verify external tools

```bash
RNAfold --version
b2ct --help
```

---

## Quick Start

### Standard SpCas9 scoring

```bash
python DGD.py input.fa
python DGD.py input.fa --output results.csv --models ./models --verbose
```

### Broad-PAM Cas9 variant scoring (xCas9, SpCas9-NG, etc.)

```bash
python Variants/DGDVar.py input.fa
python Variants/DGDVar.py input.fa --output DGDVar_results.csv --models ./models
```

### Base editor guide design (ABE / CBE)

```bash
python "Base-Editors (BE)/DGDbaseeditor.py" input.fa
python "Base-Editors (BE)/DGDbaseeditor.py" input.fa --output DGDBE_results.csv
```

### CLI arguments (all three scripts)

| Argument | Description | Default |
|---|---|---|
| `FASTA` | Input FASTA file (sequences must be 100–10,000 nt) | *(required)* |
| `--output` / `-o` | Output CSV file path | script-dependent |
| `--models` / `-m` | Directory containing trained `.h5` model files | `.` |
| `--verbose` / `-v` | Enable verbose/debug logging | off |

### Input FASTA format

```
>Gene1_sequence
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC...
>Gene2_sequence
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

Sequences must be between 100 and 10,000 nucleotides. Shorter or longer sequences are skipped with a warning.

---

## Pipeline Details (12 Steps)

All three pipeline scripts run the same 12-step workflow — the key differences are in PAM scanning strategy, guide window size, scaffold sequence, and scoring models.

```
Input FASTA
    │
    ▼ Step  1  scan_guides()
         Scan for NGG (or all-PAM) sites → Structure_file.csv
    │
    ▼ Step  2  make_fasta_for_rnafold()
         Append sgRNA scaffold to each guide → Structure_Connection.fa
    │
    ▼ Step  3  compute_target_features()
         One-hot, entropy, GC, Tm → Target_sequence_feature.csv
    │
    ▼ Step  4  RNAfold  [external]
         Predict RNA secondary structure (dot-bracket format)
    │
    ▼ Step  5  parse_rnafold_output()
         Parse b2ct connection table → Structure_out.txt
    │
    ▼ Step  6  connection_to_matrix  [C++ binary]
         Build binary base-pair matrix → Structure_basepairs.csv
    │
    ▼ Step  7  extract_spacer_scaffold_pairs()
         Filter spacer–scaffold contacts → spacer_scaffold_basepairs.csv
    │
    ▼ Step  8  compute_connection_frequency()
         Position metadata → spacer_scaffold_feature.csv
    │
    ▼ Step  9  annotate_structure_regions()
         Label positions (R, TL, AR, LR, SL1, SL2, SL3, NS) → Structural_annotation.csv
    │
    ▼ Step 10  build_features()
         Per-region connectivity + sequence features → Feature_Data_Spacer_Scaffold.csv
    │
    ▼ Step 11  compute_final_features()
         Spacer–scaffold GC ratio, Gibbs energy, monomer/dimer counts → Deep_learning_file.csv
    │
    ▼ Step 12  score_guides()
         CNN ensemble predictions → DGD.csv
```

### Structural region annotation (Step 9)

Each spacer–scaffold base-pair is assigned to one of 8 structural regions:

| Label | Region | Scaffold positions |
|---|---|---|
| R | Repeat | 21–32 |
| TL | Tetraloop | 33–36 |
| AR | Anti-repeat | 37–49 |
| LR | Lower repeat | 63–67 |
| SL1 | Stem-loop 1 | 54–58 |
| SL2 | Stem-loop 2 | 73–76 |
| SL3 | Stem-loop 3 | 88–90 |
| NS | Non-structural | all others |

---

## Variants — DGDVar.py (Broad-PAM Cas9 Variants)

`Variants/DGDVar.py` scans all 16 possible dinucleotide PAM sequences to support broad-PAM Cas9 variants such as xCas9 and SpCas9-NG. It scores each guide against 9 independent CNN ensembles:

| Output column | Cas9 variant |
|---|---|
| `DGDSpCas9` | Wild-type SpCas9 |
| `DGDeSpCas9` | eSpCas9 (enhanced specificity) |
| `DGDHypaCas9` | HypaCas9 |
| `DGDSpCas9Hf1` | SpCas9-HF1 (high-fidelity) |
| `DGDSniperCas9` | Sniper-Cas9 |
| `DGDevoCas9` | evoCas9 |
| `DGDxCas9` | xCas9 (expanded PAM) |
| `DGDSpCas9VRQR` | SpCas9-VRQR (NGA PAM) |
| `DGDSpCas9NG` | SpCas9-NG (minimal PAM) |

Scores below 0 are reported as `no_activity`. The output file is `DGDVar.csv`.

**Model directories expected under `--models`:**
```
./SpCas9/models/
./eSpCas9/models/
./HypaCas9/models/
./SpCas9-Hf1/models/
./Sniper-Cas9/models/
./evoCas9/models/
./xCas9/models/
./SpCas9-VRQR/models/
./SpCas9-NG/models/
```

---

## Base Editors — DGDbaseeditor.py

`Base-Editors (BE)/DGDbaseeditor.py` is adapted for adenine base editors (ABE) and cytosine base editors (CBE). Key differences from the standard pipeline:

| Parameter | DGD.py | DGDbaseeditor.py |
|---|---|---|
| Guide window | 30 nt | 24 nt |
| PAM | NGG only | NGG only |
| Guide extraction | positions 4:24 of window | positions 1:21 of window |
| Scaffold | SpCas9 sgRNA scaffold | ABE/CBE-optimised scaffold |
| Scaffold length | 112 nt | 102 nt |
| Models | SpCas9 DGD ensemble | ABE and CBE ensembles |
| Output columns | `DGD` | `DGDABE`, `DGDCBE` |
| Output file | `DGD.csv` | `DGDBE.csv` |

**Model directories expected under `--models`:**
```
./ABE/models/
./CBE/models/
```

---

## Accessory Scripts

Each pipeline step can also be run independently via its standalone script in the `Accessory/` folder. All scripts accept `--help` for full usage.

```bash
# Reverse complement a sequence
python Accessory/get_sequence.py --sequence ATCGATCGATCG

# Parse a FASTA file and list sequences
python Accessory/get_sequence.py --fasta input.fa

# Step 2: Build scaffold-appended FASTA
python Accessory/fastamaker.py --input Structure_file.csv --fasta-out Structure_Connection.fa

# Step 3: Compute target sequence features
python Accessory/targetsequence.py --input Structure_file.csv --output Target_sequence_feature.csv

# Step 5: Parse RNAfold output
python Accessory/Connectstructure.py --input Structure_Connection.outs --output Structure_out.txt

# Step 7: Extract spacer–scaffold base pairs
python Accessory/spacerscaffold.py --input Structure_basepairs.csv --output spacer_scaffold_basepairs.csv

# Step 8: Connection frequency table
python Accessory/spacerconnectionfrequency.py

# Step 9: Annotate structural regions
python Accessory/serialconnection.py

# Step 10: Build structural feature matrix
python Accessory/featuremaker.py

# Step 11: Final spacer–scaffold features
python Accessory/finalfeatures.py

# Step 12: Score with DGD ensemble
python Accessory/score_deep.py --models ./models --input Deep_learning_file.csv --output DGD.csv
```

---

## Output Files

### Main output — `DGD.csv`

| Column | Description |
|---|---|
| `ID` | Guide identifier: `{SeqID}:{Start}:{End}:{Strand}` |
| `Start` | Guide start position in input sequence (0-indexed) |
| `End` | Guide end position |
| `Strand` | `+` (forward) or `-` (reverse) |
| `Sequence` | 30 bp guide sequence (including PAM context) |
| `DGD` | Predicted on-target activity score (0–1) |

### Intermediate working files

| File | Generated at step |
|---|---|
| `Structure_file.csv` | Step 1 — Guide candidates |
| `Structure_Connection.fa` | Step 2 — Scaffold-appended FASTA |
| `Target_sequence_feature.csv` | Step 3 — Sequence features |
| `Structure_Connection.out` | Step 4 — RNAfold dot-bracket output |
| `Structure_out.txt` | Step 5 — Tabular connection matrix |
| `Structure_basepairs.csv` | Step 6 — Binary base-pair matrix |
| `spacer_scaffold_basepairs.csv` | Step 7 — Spacer–scaffold base pairs |
| `spacer_scaffold_feature.csv` | Step 8 — Position frequency metadata |
| `Structural_annotation.csv` | Step 9 — Region-annotated positions |
| `Feature_Data_Spacer_Scaffold.csv` | Step 10 — Combined feature matrix |
| `Deep_learning_file.csv` | Step 11 — Final model input |

---

## License

MIT

MIT License — see `LICENSE` file for details.
