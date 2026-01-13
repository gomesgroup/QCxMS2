#!/usr/bin/env python3
"""
Example: Run QCxMS2 benchmark against MassSpecGym experimental data.

This script demonstrates how to use the benchmark framework to validate
QCxMS2 calculated spectra against experimental MS/MS data from MassSpecGym.

Usage:
    python run_benchmark.py
    python run_benchmark.py --output-dir ./my_results
    python run_benchmark.py --threshold 0.7 --metrics cosine spectral_entropy

Requirements:
    pip install datasets numpy

The benchmark will:
1. Load sample molecules with simulated QCxMS2 calculated spectra
2. Search MassSpecGym (231K spectra) for matching experimental data
3. Calculate similarity metrics (cosine, spectral entropy, etc.)
4. Generate a summary report with statistics
5. Export results to JSON and CSV

Note: First run will download MassSpecGym dataset (~2GB) from Hugging Face.
"""

import argparse
import logging
import sys
from pathlib import Path

# Add parent directories to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def get_sample_molecules():
    """
    Get sample molecules with simulated QCxMS2 calculated spectra.
    
    In real usage, these would come from actual QCxMS2 calculations.
    Here we use representative EI mass spectra for demonstration.
    """
    return [
        # Ethanol (C2H6O)
        {
            "smiles": "CCO",
            "formula": "C2H6O",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "name": "Ethanol",
            "calculation_method": "EI",
            "peaks": [
                (15.0, 8.0),    # CH3+
                (26.0, 5.0),    # C2H2+
                (27.0, 15.0),   # C2H3+
                (29.0, 22.0),   # CHO+
                (31.0, 100.0),  # CH3O+ (base peak)
                (45.0, 55.0),   # C2H5O+ / CHO2+
                (46.0, 85.0),   # M+ (molecular ion)
            ],
        },
        # Acetone (C3H6O)
        {
            "smiles": "CC(=O)C",
            "formula": "C3H6O",
            "inchikey": "CSCPPACGZOOCGX-UHFFFAOYSA-N",
            "name": "Acetone",
            "calculation_method": "EI",
            "peaks": [
                (15.0, 18.0),   # CH3+
                (26.0, 8.0),    # C2H2+
                (27.0, 12.0),   # C2H3+
                (28.0, 5.0),    # CO+
                (29.0, 8.0),    # CHO+
                (42.0, 35.0),   # CH3CO+ - CO = C2H2O+
                (43.0, 100.0),  # CH3CO+ (base peak)
                (58.0, 28.0),   # M+ (molecular ion)
            ],
        },
        # Benzene (C6H6)
        {
            "smiles": "c1ccccc1",
            "formula": "C6H6",
            "inchikey": "UHOVQNZJYSORNB-UHFFFAOYSA-N",
            "name": "Benzene",
            "calculation_method": "EI",
            "peaks": [
                (26.0, 5.0),    # C2H2+
                (27.0, 7.0),    # C2H3+
                (28.0, 3.0),
                (37.0, 4.0),    # C3H+
                (38.0, 6.0),    # C3H2+
                (39.0, 15.0),   # C3H3+
                (50.0, 18.0),   # C4H2+
                (51.0, 22.0),   # C4H3+
                (52.0, 20.0),   # C4H4+
                (63.0, 5.0),    # C5H3+
                (74.0, 8.0),    # C6H2+
                (77.0, 25.0),   # C6H5+ (phenyl)
                (78.0, 100.0),  # M+ (molecular ion, base peak)
            ],
        },
        # Methanol (CH4O)
        {
            "smiles": "CO",
            "formula": "CH4O",
            "inchikey": "OKKJLVBELUTLKV-UHFFFAOYSA-N",
            "name": "Methanol",
            "calculation_method": "EI",
            "peaks": [
                (15.0, 12.0),   # CH3+
                (28.0, 8.0),    # CO+
                (29.0, 65.0),   # CHO+
                (31.0, 100.0),  # CH3O+ (base peak)
                (32.0, 72.0),   # M+ (molecular ion)
            ],
        },
        # Acetic acid (C2H4O2)
        {
            "smiles": "CC(=O)O",
            "formula": "C2H4O2",
            "inchikey": "QTBSBXVTEAMEQO-UHFFFAOYSA-N",
            "name": "Acetic acid",
            "calculation_method": "EI",
            "peaks": [
                (15.0, 20.0),   # CH3+
                (28.0, 8.0),    # CO+
                (29.0, 15.0),   # CHO+
                (42.0, 12.0),   # CH2CO+
                (43.0, 100.0),  # CH3CO+ (base peak)
                (45.0, 85.0),   # CHO2+
                (60.0, 55.0),   # M+ (molecular ion)
            ],
        },
        # Toluene (C7H8)
        {
            "smiles": "Cc1ccccc1",
            "formula": "C7H8",
            "inchikey": "YXFVVABEGXRONW-UHFFFAOYSA-N",
            "name": "Toluene",
            "calculation_method": "EI",
            "peaks": [
                (27.0, 8.0),
                (39.0, 18.0),   # C3H3+
                (50.0, 12.0),   # C4H2+
                (51.0, 15.0),   # C4H3+
                (63.0, 10.0),   # C5H3+
                (65.0, 18.0),   # C5H5+
                (77.0, 8.0),    # C6H5+
                (89.0, 5.0),    # C7H5+
                (91.0, 100.0),  # C7H7+ (tropylium, base peak)
                (92.0, 70.0),   # M+ (molecular ion)
            ],
        },
        # Phenol (C6H6O)
        {
            "smiles": "Oc1ccccc1",
            "formula": "C6H6O",
            "inchikey": "ISWSIDIOOBJBQZ-UHFFFAOYSA-N",
            "name": "Phenol",
            "calculation_method": "EI",
            "peaks": [
                (28.0, 5.0),    # CO+
                (38.0, 8.0),
                (39.0, 20.0),   # C3H3+
                (40.0, 8.0),
                (50.0, 10.0),
                (51.0, 8.0),
                (55.0, 12.0),
                (63.0, 5.0),
                (65.0, 20.0),   # C5H5+
                (66.0, 25.0),   # C5H6+ / C4H2O+
                (94.0, 100.0),  # M+ (molecular ion, base peak)
            ],
        },
        # Cyclohexane (C6H12)
        {
            "smiles": "C1CCCCC1",
            "formula": "C6H12",
            "inchikey": "XDTMQSROBMDMFD-UHFFFAOYSA-N",
            "name": "Cyclohexane",
            "calculation_method": "EI",
            "peaks": [
                (27.0, 25.0),
                (28.0, 10.0),
                (29.0, 15.0),
                (39.0, 30.0),
                (41.0, 55.0),
                (42.0, 25.0),
                (43.0, 15.0),
                (55.0, 40.0),
                (56.0, 100.0),  # Base peak (M-28, loss of C2H4)
                (69.0, 25.0),
                (84.0, 70.0),   # M+ (molecular ion)
            ],
        },
        # Aniline (C6H7N)
        {
            "smiles": "Nc1ccccc1",
            "formula": "C6H7N",
            "inchikey": "PAYRUJLWNCNPSJ-UHFFFAOYSA-N",
            "name": "Aniline",
            "calculation_method": "EI",
            "peaks": [
                (27.0, 8.0),
                (28.0, 5.0),
                (38.0, 5.0),
                (39.0, 18.0),
                (40.0, 8.0),
                (50.0, 10.0),
                (51.0, 12.0),
                (52.0, 8.0),
                (63.0, 8.0),
                (64.0, 8.0),
                (65.0, 18.0),
                (66.0, 30.0),
                (93.0, 100.0),  # M+ (molecular ion, base peak)
            ],
        },
        # Formic acid (CH2O2)
        {
            "smiles": "O=CO",
            "formula": "CH2O2",
            "inchikey": "BDAGIHXWWSANSR-UHFFFAOYSA-N",
            "name": "Formic acid",
            "calculation_method": "EI",
            "peaks": [
                (17.0, 8.0),    # OH+
                (28.0, 15.0),   # CO+
                (29.0, 100.0),  # CHO+ (base peak)
                (44.0, 5.0),    # CO2+
                (45.0, 45.0),   # CHO2+
                (46.0, 70.0),   # M+ (molecular ion)
            ],
        },
    ]


def progress_callback(current: int, total: int, smiles: str):
    """Progress callback for benchmark."""
    print(f"\rBenchmarking: {current}/{total} - {smiles[:40]:<40}", end="", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Run QCxMS2 benchmark against MassSpecGym experimental data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./benchmark_results",
        help="Directory for output files (default: ./benchmark_results)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Similarity threshold for matching (default: 0.5)",
    )
    parser.add_argument(
        "--mz-tolerance",
        type=float,
        default=0.5,
        help="m/z tolerance for peak matching in Da (default: 0.5)",
    )
    parser.add_argument(
        "--metrics",
        nargs="+",
        default=["cosine", "dot_product", "spectral_entropy"],
        help="Similarity metrics to calculate (default: cosine dot_product spectral_entropy)",
    )
    parser.add_argument(
        "--max-matches",
        type=int,
        default=50,
        help="Maximum experimental matches per molecule (default: 50)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress output",
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    print("=" * 70)
    print("QCxMS2 Benchmark Suite - MassSpecGym Validation")
    print("=" * 70)
    print()
    
    # Import benchmark module
    try:
        from benchmark import QCxMS2Benchmark, BenchmarkConfig
    except ImportError as e:
        logger.error(f"Failed to import benchmark module: {e}")
        logger.error("Make sure you're running from the database_integration directory")
        sys.exit(1)
    
    # Create configuration
    config = BenchmarkConfig(
        similarity_threshold=args.threshold,
        mz_tolerance=args.mz_tolerance,
        metrics=args.metrics,
        max_experimental_matches=args.max_matches,
        output_dir=args.output_dir,
    )
    
    print("Configuration:")
    print(f"  - Similarity threshold: {config.similarity_threshold}")
    print(f"  - m/z tolerance: {config.mz_tolerance} Da")
    print(f"  - Metrics: {', '.join(config.metrics)}")
    print(f"  - Max experimental matches: {config.max_experimental_matches}")
    print(f"  - Output directory: {config.output_dir}")
    print()
    
    # Create benchmark
    benchmark = QCxMS2Benchmark(config)
    
    # Add sample molecules
    print("Loading sample molecules...")
    molecules = get_sample_molecules()
    
    for mol in molecules:
        benchmark.add_molecule(
            smiles=mol["smiles"],
            formula=mol["formula"],
            inchikey=mol["inchikey"],
            calculated_spectrum=mol["peaks"],
            name=mol["name"],
            calculation_method=mol["calculation_method"],
        )
    
    print(f"Added {benchmark.num_molecules} molecules to benchmark")
    print()
    
    # Run benchmark
    print("Running benchmark...")
    print("(First run downloads MassSpecGym dataset ~2GB from Hugging Face)")
    print()
    
    callback = None if args.quiet else progress_callback
    
    try:
        results = benchmark.run_benchmark(progress_callback=callback)
    except Exception as e:
        logger.error(f"Benchmark failed: {e}")
        print()
        print("If you see an ImportError for 'datasets', install it with:")
        print("  pip install datasets")
        sys.exit(1)
    
    print()  # New line after progress
    print()
    
    # Print summary
    print("=" * 70)
    print("BENCHMARK RESULTS SUMMARY")
    print("=" * 70)
    print()
    print(f"Total molecules benchmarked: {results.total_molecules}")
    print(f"Molecules with experimental matches: {results.molecules_with_matches}")
    print(f"Molecules above similarity threshold: {results.molecules_above_threshold}")
    print(f"Total benchmark time: {results.total_time:.2f} seconds")
    print()
    
    if results.cosine_stats:
        stats = results.cosine_stats
        print("Cosine Similarity Statistics (Best Match per Molecule):")
        print(f"  Mean:   {stats.mean:.4f} +/- {stats.std:.4f}")
        print(f"  Median: {stats.median:.4f}")
        print(f"  Range:  [{stats.min_value:.4f}, {stats.max_value:.4f}]")
        print(f"  Above 0.7: {stats.above_threshold_70:.1f}%")
        print(f"  Above 0.9: {stats.above_threshold_90:.1f}%")
        print()
    
    # Print individual results
    print("Individual Molecule Results:")
    print("-" * 70)
    print(f"{'Molecule':<20} {'Formula':<10} {'Best Cosine':<12} {'Matches':<10} {'Quality':<10}")
    print("-" * 70)
    
    for mol in sorted(results.molecules, key=lambda m: m.best_cosine, reverse=True):
        quality = mol.best_match.quality_tier if mol.best_match else "N/A"
        print(f"{mol.name or mol.smiles[:20]:<20} {mol.formula:<10} {mol.best_cosine:<12.4f} {mol.num_experimental_found:<10} {quality:<10}")
    
    print()
    
    # Generate report
    print("Generating report...")
    report_path = Path(args.output_dir) / "benchmark_report.md"
    benchmark.generate_report(report_path)
    print(f"Report saved to: {report_path}")
    
    # Export results
    json_path = Path(args.output_dir) / "benchmark_results.json"
    csv_path = Path(args.output_dir) / "benchmark_results.csv"
    
    benchmark.export_results(json_path)
    benchmark.export_results(csv_path)
    
    print(f"JSON results saved to: {json_path}")
    print(f"CSV results saved to: {csv_path}")
    print()
    
    print("=" * 70)
    print("Benchmark complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
