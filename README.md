# DGD-Cas9 — Deep Guide Designer for CRISPR-Cas9

This branch (`main`) is a lightweight landing page: it shows the
**Snakemake workflow** and a small **demo input** so you can see how the
pipeline is organized, without cloning the full source tree, trained
models, or PDB-scale data.

**This branch is not runnable on its own** — the Snakefile here calls into
Python pipeline scripts (`DGD.py`, `Accessory/*.py`, etc.) that
intentionally live on the `full-pipeline` branch instead, to keep `main`
small and easy to browse.

## To actually run the pipeline

```bash
git clone https://github.com/GuideDesigner/DGD-Cas9.git
cd DGD-Cas9
git checkout full-pipeline
pip install -r requirements.txt
pip install snakemake "pulp<2.8"
snakemake --cores 1 -n     # preview the workflow
snakemake --cores 1        # run it against demo/input.fa
```

See the `full-pipeline` branch's README for full installation instructions,
the 12-step pipeline description, the three pipeline variants (standard
SpCas9, broad-PAM Cas9 variants, base editors), known issues, and CI/CD
details.

## What's here

```
main/
├── Snakefile              # The real workflow definition (same file as full-pipeline)
├── config/config.yaml     # Sample/variant configuration
├── envs/environment.yaml  # Conda environment for --use-conda
└── demo/input.fa          # Small synthetic FASTA for trying the workflow
```

## Overview

**DGD-Cas9** is a deep-learning pipeline for designing and scoring
CRISPR-Cas9 sgRNA guide candidates. It integrates sequence features, RNA
secondary structure, and spacer–scaffold base-pairing connectivity into a
CNN ensemble model to predict guide activity.

**Authors:** Vipin Menon, Jang-il Sohn, Seokju Park, Jin-Wu Nam
**Lab:** Bioinformatics & Genomics Lab (BIG Lab), Hanyang University, Seoul 04763, Korea
**Contact:** a.vipin.menon@gmail.com | jwnam@hanyang.ac.kr

## License

MIT License — see the `full-pipeline` branch's `LICENSE` file for details.
