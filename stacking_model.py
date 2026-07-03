#!/usr/bin/env python3
"""
stacking_model.py — RNA dinucleotide stacking free-energy model for DGD
========================================================================
Author : Vipin Menon, BIG Lab, Hanyang University
Date   : August 2021 (modernized 2026)

Provides `return_stacking_model()` which returns a dictionary-of-dictionaries
mapping RNA dinucleotide pairs to their stacking free energies (kcal/mol).
"""


def return_stacking_model() -> dict:
    """
    Return the RNA stacking free-energy model.

    Constructs and returns a nested dictionary mapping dinucleotide pairs
    to their stacking free energies in kcal/mol at 37°C.

    ****************************************************************
    KNOWN ISSUE -- INCOMPLETE TABLE, NEEDS ATTENTION BEFORE TRUSTING
    RESULTS THAT DEPEND ON Spacer_Scaffold_Gibbs_Energy:
    ****************************************************************
    A standard nearest-neighbour RNA stacking table has 16 top-level
    dinucleotide contexts (AA, AU, AG, AC, UA, UU, UG, UC, GA, GU, GG,
    GC, CA, CU, CG, CC). This dict only defines 4 of them (AA, AU, AG,
    AC). Since real guide/scaffold sequences routinely contain the other
    12 dinucleotides, DGD.py's _calculate_stacking_energy() used to raise
    a hard KeyError the first time it looked one up -- i.e. step 11
    (compute_final_features) could not complete for any real input.

    Rather than invent the missing thermodynamic values (they should come
    from the original published parameter set this model was trained
    against, not be guessed), calculate_stacking_energy() below now falls
    back to 0.0 kcal/mol for any missing context and logs a warning. This
    keeps the pipeline runnable end-to-end for development/testing, but
    the resulting Spacer_Scaffold_Gibbs_Energy feature is only fully
    correct for guides whose paired dimers happen to only use AA/AU/AG/AC
    contexts. Please supply the complete 16-context table before relying
    on this feature for real predictions.

    Returns:
        A dict where keys are dinucleotide strings (e.g., 'AA', 'GC')
        and values are dicts mapping base-pair contexts to stacking energies.
    """
    stacking_model = {
        'AA': {
            'AA': -0.93, 'AU': -0.88, 'AG': -0.58, 'AC': -0.52,
            'UA': -0.93, 'UU': -0.80, 'UG': -0.50, 'UC': -0.46,
            'GA': -0.59, 'GU': -0.50, 'GG': -0.29, 'GC': -0.12,
            'CA': -0.56, 'CU': -0.48, 'CG': -0.32, 'CC': -0.29
        },
        'AU': {
            'AA': -0.95, 'AU': -0.84, 'AG': -0.64, 'AC': -0.52,
            'UA': -0.95, 'UU': -0.79, 'UG': -0.61, 'UC': -0.47,
            'GA': -0.65, 'GU': -0.60, 'GG': -0.47, 'GC': -0.32,
            'CA': -0.62, 'CU': -0.55, 'CG': -0.48, 'CC': -0.35
        },
        'AG': {
            'AA': -0.89, 'AU': -0.84, 'AG': -0.56, 'AC': -0.51,
            'UA': -0.89, 'UU': -0.78, 'UG': -0.48, 'UC': -0.44,
            'GA': -0.57, 'GU': -0.48, 'GG': -0.27, 'GC': -0.10,
            'CA': -0.54, 'CU': -0.46, 'CG': -0.30, 'CC': -0.27
        },
        'AC': {
            'AA': -0.89, 'AU': -0.80, 'AG': -0.54, 'AC': -0.44,
            'UA': -0.89, 'UU': -0.73, 'UG': -0.44, 'UC': -0.35,
            'GA': -0.55, 'GU': -0.44, 'GG': -0.22, 'GC': -0.04,
            'CA': -0.52, 'CU': -0.42, 'CG': -0.24, 'CC': -0.14
        }
    }
    return stacking_model
