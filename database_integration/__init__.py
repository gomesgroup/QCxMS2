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
    
    # Benchmark against MassSpecGym (231K spectra)
    from qcxms2.database_integration.benchmark import QCxMS2Benchmark, BenchmarkConfig
    benchmark = QCxMS2Benchmark(BenchmarkConfig(similarity_threshold=0.7))
    benchmark.add_molecule(smiles="CCO", formula="C2H6O", inchikey="...", calculated_spectrum=[...])
    results = benchmark.run_benchmark()
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

# QCxMS2 output parsers
from .parsers import (
    QCxMS2Result,
    QCxMS2FragmentInfo,
    parse_qcxms2_results,
    parse_qcxms2_spectrum,
    parse_qcxms2_allpeaks,
    parse_qcxms2_output_dir,
)

# Benchmark framework
from .benchmark import (
    QCxMS2Benchmark,
    BenchmarkConfig,
    BenchmarkResult,
    MoleculeBenchmark,
)

__version__ = "1.2.0"

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
    # Benchmark
    "QCxMS2Benchmark",
    "BenchmarkConfig",
    "BenchmarkResult",
    "MoleculeBenchmark",
]
