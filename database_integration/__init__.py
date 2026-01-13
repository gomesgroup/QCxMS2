"""
QCxMS2 Database Integration Module

This module provides connectors to experimental mass spectrometry databases,
enabling QCxMS2 to compare calculated spectra with experimental data.

Architecture inspired by Explorer's database_integration system.

Supported Databases:
- MassBank (massbank.eu): 40K+ mass spectra with REST API
- MassSpecGym: NeurIPS 2024 benchmark with 231K curated MS/MS spectra
- SDBS (sdbs.db.aist.go.jp): IR, NMR, MS, Raman for organic compounds
- NIST MS (webbook.nist.gov): NIST Chemistry WebBook mass spectra
- MoNA (mona.fiehnlab.ucdavis.edu): MassBank of North America

Usage:
    from qcxms2.database_integration import create_ms_connector, DatabaseType
    
    # Create connector
    massbank = create_ms_connector(DatabaseType.MASSBANK)
    
    # Or use MassSpecGym (231K curated spectra, downloaded via Hugging Face)
    massspecgym = create_ms_connector(DatabaseType.MASSSPECGYM)
    
    # Search for spectrum by formula
    from qcxms2.database_integration import SpectrumSearchCriteria
    criteria = SpectrumSearchCriteria(formula="C6H12O6", ionization_mode="EI")
    results = massbank.search_spectrum(criteria)
    
    # Compare with experimental data
    from qcxms2.database_integration import compare_spectra
    similarity = compare_spectra(calculated_peaks, experimental_peaks)
"""

from .base import MSDatabaseConnector, ExperimentalSpectrum
from .config import DatabaseConfig
from .connectors import (
    MassBankConnector,
    MassSpecGymConnector,
    MoNAConnector,
    NISTMSConnector,
    SDBSConnector,
)
from .data_classes import (
    ExperimentalSpectrum,
    MSSearchResult,
    Peak,
    SpectrumComparisonResult,
)
from .enums import DatabaseType, IonizationMode, SpectrumType
from .factory import create_ms_connector, create_massspecgym_connector
from .search_criteria import SpectrumSearchCriteria
from .utils import compare_spectra, parse_spectrum_file

__version__ = "1.1.0"

__all__ = [
    # Base classes
    "MSDatabaseConnector",
    "ExperimentalSpectrum",
    # Config
    "DatabaseConfig",
    # Connectors
    "MassBankConnector",
    "MassSpecGymConnector",
    "SDBSConnector",
    "NISTMSConnector",
    "MoNAConnector",
    # Data classes
    "ExperimentalSpectrum",
    "MSSearchResult",
    "Peak",
    "SpectrumComparisonResult",
    # Enums
    "DatabaseType",
    "IonizationMode",
    "SpectrumType",
    # Factory
    "create_ms_connector",
    "create_massspecgym_connector",
    # Search
    "SpectrumSearchCriteria",
    # Utils
    "compare_spectra",
    "parse_spectrum_file",
]
