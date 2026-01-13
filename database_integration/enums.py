"""
Enumerations for mass spectrometry database integration.

This module defines enum types for MS databases, ionization modes, and spectrum types.
"""

from enum import Enum


class DatabaseType(Enum):
    """Supported mass spectrometry databases."""

    MASSBANK = "massbank"  # MassBank EU - massbank.eu
    MONA = "mona"  # MassBank of North America - mona.fiehnlab.ucdavis.edu
    NIST_MS = "nist_ms"  # NIST Chemistry WebBook - webbook.nist.gov
    SDBS = "sdbs"  # Spectral Database for Organic Compounds - sdbs.db.aist.go.jp
    GNPS = "gnps"  # Global Natural Products Social Molecular Networking
    METLIN = "metlin"  # METLIN Metabolomics Database
    HMDB = "hmdb"  # Human Metabolome Database
    LIPIDMAPS = "lipidmaps"  # LIPID MAPS Structure Database
    CUSTOM = "custom"  # Custom/local database


class IonizationMode(Enum):
    """Ionization modes for mass spectrometry."""

    EI = "EI"  # Electron Ionization (Electron Impact)
    CID = "CID"  # Collision-Induced Dissociation
    DEA = "DEA"  # Dissociative Electron Attachment
    ESI_POSITIVE = "ESI+"  # Electrospray Ionization (positive)
    ESI_NEGATIVE = "ESI-"  # Electrospray Ionization (negative)
    APCI_POSITIVE = "APCI+"  # Atmospheric Pressure Chemical Ionization (positive)
    APCI_NEGATIVE = "APCI-"  # Atmospheric Pressure Chemical Ionization (negative)
    MALDI = "MALDI"  # Matrix-Assisted Laser Desorption/Ionization
    CI = "CI"  # Chemical Ionization
    FI = "FI"  # Field Ionization
    FD = "FD"  # Field Desorption
    FAB = "FAB"  # Fast Atom Bombardment
    OTHER = "Other"


class SpectrumType(Enum):
    """Type of mass spectrum."""

    MS1 = "MS1"  # Single MS (parent ions)
    MS2 = "MS2"  # MS/MS (fragment ions)
    MSN = "MSn"  # Multiple stages of MS
    THEORETICAL = "Theoretical"  # Calculated/predicted spectrum
    EXPERIMENTAL = "Experimental"  # Measured spectrum


class SearchType(Enum):
    """Type of spectrum search/matching."""

    IDENTITY = "identity"  # Exact match by InChI or structure
    SIMILARITY = "similarity"  # Similarity search by spectrum
    SUBSTRUCTURE = "substructure"  # Search by substructure
    FORMULA = "formula"  # Search by molecular formula
    EXACT_MASS = "exact_mass"  # Search by exact mass
    PRECURSOR_MZ = "precursor_mz"  # Search by precursor m/z


class SimilarityMetric(Enum):
    """Metrics for spectrum similarity comparison."""

    COSINE = "cosine"  # Cosine similarity
    DOT_PRODUCT = "dot_product"  # Dot product
    REVERSE_DOT_PRODUCT = "reverse_dot_product"  # Reverse dot product
    EUCLIDEAN = "euclidean"  # Euclidean distance
    MANHATTAN = "manhattan"  # Manhattan distance
    SPECTRAL_ENTROPY = "spectral_entropy"  # Spectral entropy similarity
    MODIFIED_COSINE = "modified_cosine"  # Modified cosine for shifted peaks
