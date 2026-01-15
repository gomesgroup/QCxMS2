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

__all__ = [
    "FragmentPrioritizer",
    "PrioritizedFragment",
    "prioritize_fragments_from_spectrum",
    "prioritize_fragments_from_smiles",
    "EarlyStoppingMonitor",
    "ConvergenceStatus",
    "monitor_calculation",
]
