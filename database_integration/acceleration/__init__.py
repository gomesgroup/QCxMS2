"""
QCxMS2 Acceleration Module

Provides ML-based acceleration for QCxMS2 calculations using:
- DreaMS embeddings for fragment prioritization
- MEDUSA LSTM for formula prediction
- Cached barriers for repeated calculations

Key Components:
- FragmentPrioritizer: Use DreaMS to predict expected fragments
- FormulaPruner: Use MEDUSA to constrain fragmentation search (planned)
- EarlyStopping: Stop calculation when spectrum converges (planned)
"""

from .fragment_prioritizer import (
    FragmentPrioritizer,
    PrioritizedFragment,
    prioritize_fragments_from_spectrum,
    prioritize_fragments_from_smiles,
)

__all__ = [
    "FragmentPrioritizer",
    "PrioritizedFragment",
    "prioritize_fragments_from_spectrum",
    "prioritize_fragments_from_smiles",
]
