"""
DGD-Cas9 Snakemake workflow
=============================
Wraps the 12-step DGD guide-scoring pipeline. For the "standard" SpCas9
variant, each step is wired as its own Snakemake rule (using the
Accessory/*.py standalone scripts) so you can see -- and resume from --
exactly where the pipeline is, and so it can run as far as steps 1-11
(guide scanning through feature engineering) without needing trained
model weights. The "broad_pam" and "base_editor" variants are wired as a
single rule each, matching their all-in-one CLI usage.

Two things this repo does NOT include (see README "Known Issues"):
  - connection_to_matrix.cpp (C++ source for Step 6) -- build it with
    `make` once you have the source; until then, Step 6 fails clearly.
  - Trained .h5 model files (models/) -- needed only for the final
    scoring step, controlled by config['run_scoring'].

Usage:
    snakemake --cores 1 -n              # preview the DAG
    snakemake --cores 1                 # run with deps already on PATH
    snakemake --cores 1 --use-conda     # run using envs/environment.yaml
"""

configfile: "config/config.yaml"

FASTA        = config["fasta"]
VARIANT      = config.get("variant", "standard")
MODELS_DIR   = config.get("models_dir", "models")
RUN_SCORING  = config.get("run_scoring", False)

VALID_VARIANTS = ("standard", "broad_pam", "base_editor")
if VARIANT not in VALID_VARIANTS:
    raise ValueError(f"config['variant'] must be one of {VALID_VARIANTS}, got '{VARIANT}'")


def _standard_targets():
    targets = ["Deep_learning_file.csv"]
    if RUN_SCORING:
        targets.append("DGD.csv")
    return targets


rule all:
    input:
        _standard_targets() if VARIANT == "standard"
        else ["DGDVar.csv"] if VARIANT == "broad_pam"
        else ["DGDBE.csv"],


# ===========================================================================
# Standard SpCas9 variant (DGD.py) -- granular, step-by-step rules
# ===========================================================================

rule scan_guides:
    """Step 1: scan the input FASTA for NGG-PAM guide candidates."""
    input:
        fasta=FASTA,
    output:
        "Structure_file.csv",
    log:
        "logs/01_scan_guides.log",
    conda:
        "envs/environment.yaml"
    shell:
        "mkdir -p logs && "
        "python -c \"from DGD import scan_guides; scan_guides('{input.fasta}')\" > {log} 2>&1"


rule build_rnafold_fasta:
    """Step 2: append the sgRNA scaffold to each guide for RNAfold."""
    input:
        "Structure_file.csv",
    output:
        fasta="Structure_Connection.fa",
        csv="Structure_Connection.csv",
    log:
        "logs/02_fastamaker.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/fastamaker.py --input {input} "
        "--fasta-out {output.fasta} --csv-out {output.csv} > {log} 2>&1"


rule target_features:
    """Step 3: sequence/thermodynamic/compositional features per guide."""
    input:
        "Structure_file.csv",
    output:
        "Target_sequence_feature.csv",
    log:
        "logs/03_targetsequence.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/targetsequence.py --input {input} --output {output} > {log} 2>&1"


rule rnafold:
    """Step 4 (external tools): RNAfold + b2ct secondary-structure prediction."""
    input:
        "Structure_Connection.fa",
    output:
        out="Structure_Connection.out",
        outs="Structure_Connection.outs",
    log:
        "logs/04_rnafold.log",
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p logs
        (RNAfold -j0 --noPS < {input} > {output.out}) > {log} 2>&1
        (b2ct < {output.out} > {output.outs}) >> {log} 2>&1
        """


rule parse_rnafold:
    """Step 5: parse b2ct output into a tabular connection matrix."""
    input:
        "Structure_Connection.outs",
    output:
        "Structure_out.txt",
    log:
        "logs/05_connectstructure.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/Connectstructure.py --input {input} --output {output} > {log} 2>&1"


rule connection_to_matrix:
    """
    Step 6 (external C++ binary): build the full base-pair connection matrix.

    connection_to_matrix.cpp is not included in this repository (see
    README "Known Issues") -- this rule will fail with a clear message
    until you supply the source and run `make`.
    """
    input:
        "Structure_out.txt",
    output:
        "Structure_basepairs.csv",
    log:
        "logs/06_connection_to_matrix.log",
    shell:
        """
        mkdir -p logs
        if [ ! -x ./connection_to_matrix ]; then
            echo "ERROR: ./connection_to_matrix binary not found or not executable." | tee {log}
            echo "This is built from connection_to_matrix.cpp via 'make' -- that" | tee -a {log}
            echo "source file is not present in this branch's history (see README" | tee -a {log}
            echo "'Known Issues'). Supply it and run 'make' before this rule can run." | tee -a {log}
            exit 1
        fi
        ./connection_to_matrix {input} 102 > {output} 2> {log}
        """


rule spacer_scaffold_pairs:
    """Step 7: extract spacer-to-scaffold base-pair columns."""
    input:
        "Structure_basepairs.csv",
    output:
        "spacer_scaffold_basepairs.csv",
    log:
        "logs/07_spacerscaffold.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/spacerscaffold.py --input {input} --output {output} > {log} 2>&1"


rule connection_frequency:
    """Step 8: pivot spacer-scaffold connections, extract position labels."""
    input:
        "spacer_scaffold_basepairs.csv",
    output:
        "spacer_scaffold_feature.csv",
    log:
        "logs/08_spacerconnectionfrequency.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/spacerconnectionfrequency.py --input {input} --output {output} > {log} 2>&1"


rule annotate_regions:
    """Step 9: label each connection with its structural region (R/TL/AR/...)."""
    input:
        "spacer_scaffold_feature.csv",
    output:
        "Structural_annotation.csv",
    log:
        "logs/09_serialconnection.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/serialconnection.py --input {input} --output {output} > {log} 2>&1"


rule build_features:
    """Step 10: combine sequence features with structural connectivity features."""
    input:
        bp="Structure_basepairs.csv",
        seq="Target_sequence_feature.csv",
        annot="Structural_annotation.csv",
    output:
        "Feature_Data_Spacer_Scaffold.csv",
    log:
        "logs/10_featuremaker.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/featuremaker.py --bp-input {input.bp} --seq-input {input.seq} "
        "--annot-input {input.annot} --output {output} > {log} 2>&1"


rule final_features:
    """Step 11: monomer/dimer counts + stacking energy -> final model input."""
    input:
        rnafold_out="Structure_Connection.out",
        features="Feature_Data_Spacer_Scaffold.csv",
    output:
        "Deep_learning_file.csv",
    log:
        "logs/11_finalfeatures.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/finalfeatures.py --rnafold-out {input.rnafold_out} "
        "--feature-input {input.features} --output {output} > {log} 2>&1"


rule score_standard:
    """Step 12: score guides with the trained DGD CNN ensemble."""
    input:
        features="Deep_learning_file.csv",
        struct="Structure_file.csv",
    output:
        "DGD.csv",
    log:
        "logs/12_score_deep.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Accessory/score_deep.py --models {MODELS_DIR} --input {input.features} "
        "--struct-input {input.struct} --output {output} > {log} 2>&1"


# ===========================================================================
# Broad-PAM and Base-Editor variants -- single-rule, full-CLI invocation
# (these scripts don't have per-step Accessory wrappers, so we call the
# documented end-to-end CLI directly, same as a user would by hand)
# ===========================================================================

rule score_broad_pam:
    """Variants/DGDVar.py: full pipeline for 9 broad-PAM Cas9 variants."""
    input:
        fasta=FASTA,
    output:
        "DGDVar.csv",
    log:
        "logs/dgdvar_full_run.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Variants/DGDVar.py {input.fasta} --output {output} "
        "--models {MODELS_DIR} > {log} 2>&1"


rule score_base_editor:
    """Base-Editors/DGDbaseeditor.py: full pipeline for ABE/CBE base editors."""
    input:
        fasta=FASTA,
    output:
        "DGDBE.csv",
    log:
        "logs/dgdbe_full_run.log",
    conda:
        "envs/environment.yaml"
    shell:
        "python Base-Editors/DGDbaseeditor.py {input.fasta} --output {output} "
        "--models {MODELS_DIR} > {log} 2>&1"
