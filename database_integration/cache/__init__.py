"""
Fragment Cache Module for QCxMS2

Provides caching of computed fragments and their NEB barriers to accelerate
repeated calculations on similar molecules.

Key Components:
- FragmentCache: Main cache interface using SQLite + m/z indexing
- FragmentEntry: Data class for cached fragments
- QCxMS2OutputParser: Parser for extracting fragments from QCxMS2 output

Usage:
    from database_integration.cache import FragmentCache
    
    cache = FragmentCache("/path/to/cache.db")
    cache.seed_from_qcxms2_output("/path/to/calculation")
    
    # Later: lookup cached barrier
    result = cache.find_fragment(mz=77.039, formula="C6H5")
    if result:
        print(f"Cached barrier: {result.barrier_kcal_mol} kcal/mol")
"""

from .fragment_cache import FragmentCache, FragmentEntry
from .qcxms2_output_parser import QCxMS2OutputParser

__all__ = ["FragmentCache", "FragmentEntry", "QCxMS2OutputParser"]
