"""
QCxMS2 Acceleration Module

Provides ML-based acceleration for QCxMS2 calculations using:
- DreaMS embeddings for fragment prioritization
- DreaMS-based early stopping when spectrum converges
- MEDUSA LSTM for formula prediction (planned)
- Cached barriers for repeated calculations

Key Components:
- FragmentPrioritizer: Use DreaMS to predict expected fragments
- EarlyStoppingMonitor: Stop calculation when spectrum converges
- FormulaPruner: Use MEDUSA to constrain fragmentation search (planned)
"""

from .fragment_prioritizer import (
    FragmentPrioritizer,
    PrioritizedFragment,
    prioritize_fragments_from_spectrum,
    prioritize_fragments_from_smiles,
)

from .early_stopping import (
    EarlyStoppingMonitor,
    ConvergenceStatus,
    monitor_calculation,
)

from .formula_pruning import (
    FormulaPruner,
    FormulaCandidate,
    enumerate_formulas,
    parse_formula,
    formula_to_string,
    calculate_mass,
    calculate_dbe,
    validate_valence,
    validate_nitrogen_rule,
    is_subformula,
)

__all__ = [
    # Fragment prioritization
    "FragmentPrioritizer",
    "PrioritizedFragment",
    "prioritize_fragments_from_spectrum",
    "prioritize_fragments_from_smiles",
    # Early stopping
    "EarlyStoppingMonitor",
    "ConvergenceStatus",
    "monitor_calculation",
    # Formula pruning
    "FormulaPruner",
    "FormulaCandidate",
    "enumerate_formulas",
    "parse_formula",
    "formula_to_string",
    "calculate_mass",
    "calculate_dbe",
    "validate_valence",
    "validate_nitrogen_rule",
    "is_subformula",
]
