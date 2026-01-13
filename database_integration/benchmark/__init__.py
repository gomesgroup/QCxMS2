"""
QCxMS2 Benchmark Module

A comprehensive benchmark framework for validating QCxMS2 calculated spectra
against MassSpecGym experimental data (231K MS/MS spectra from 29K molecules).

Features:
- Batch processing of multiple molecules
- Multiple similarity metrics (cosine, dot product, spectral entropy)
- Statistical analysis and quality assessment
- Export results in CSV, JSON, and markdown formats
- Summary statistics for publication

Usage:
    from qcxms2.database_integration.benchmark import (
        QCxMS2Benchmark,
        BenchmarkConfig,
        BenchmarkResult,
    )
    
    # Create benchmark with custom configuration
    config = BenchmarkConfig(
        similarity_threshold=0.7,
        metrics=["cosine", "spectral_entropy"],
        output_dir="./benchmark_results",
    )
    benchmark = QCxMS2Benchmark(config)
    
    # Add molecules to benchmark
    benchmark.add_molecule(
        smiles="CCO",
        formula="C2H6O",
        inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        calculated_spectrum={"peaks": [(15.0, 10.0), (31.0, 100.0), (45.0, 50.0)]},
    )
    
    # Run benchmark against MassSpecGym
    results = benchmark.run_benchmark()
    
    # Generate report and export
    benchmark.generate_report()
    benchmark.export_results("benchmark_results.json")
"""

from .benchmark_suite import (
    BenchmarkConfig,
    BenchmarkResult,
    MoleculeBenchmark,
    QCxMS2Benchmark,
)
from .metrics import (
    MetricsCalculator,
    Peak,
    QCxMS2ValidationMetrics,
    SimilarityMetrics,
    StatisticalAnalysis,
    calculate_cosine_similarity,
    calculate_dot_product,
    calculate_fragmentation_pattern_match,
    calculate_intensity_rank_correlation,
    calculate_mass_accuracy_score,
    calculate_qcxms2_validation_metrics,
    calculate_spectral_entropy_similarity,
    calculate_weighted_cosine,
)

__all__ = [
    # Main classes
    "QCxMS2Benchmark",
    "BenchmarkConfig",
    "BenchmarkResult",
    "MoleculeBenchmark",
    # Metrics
    "MetricsCalculator",
    "Peak",
    "QCxMS2ValidationMetrics",
    "SimilarityMetrics",
    "StatisticalAnalysis",
    # Functions
    "calculate_cosine_similarity",
    "calculate_dot_product",
    "calculate_fragmentation_pattern_match",
    "calculate_intensity_rank_correlation",
    "calculate_mass_accuracy_score",
    "calculate_qcxms2_validation_metrics",
    "calculate_spectral_entropy_similarity",
    "calculate_weighted_cosine",
]
