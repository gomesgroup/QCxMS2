"""
QCxMS2 Benchmark Suite for validating calculated spectra against MassSpecGym.

This module provides the main benchmark framework for systematically comparing
QCxMS2 calculated mass spectra against experimental data from MassSpecGym
(231K curated MS/MS spectra from 29K molecules).

Features:
- Batch processing of multiple molecules
- Multiple similarity metrics for comprehensive comparison
- Automatic matching to experimental spectra by InChI key, formula, or SMILES
- Statistical analysis and quality assessment
- Export results in JSON, CSV, and markdown formats
"""

import csv
import json
import logging
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from .metrics import (
    MetricsCalculator,
    Peak,
    SimilarityMetrics,
    StatisticalAnalysis,
    assess_match_quality,
    compute_statistical_analysis,
)

logger = logging.getLogger(__name__)


@dataclass
class BenchmarkConfig:
    """
    Configuration for QCxMS2 benchmark.
    
    Attributes:
        similarity_threshold: Minimum similarity score to consider a match (0-1)
        mz_tolerance: m/z tolerance for peak matching (Da)
        metrics: List of metrics to calculate ("cosine", "dot_product", etc.)
        max_experimental_matches: Maximum experimental spectra to compare per molecule
        output_dir: Directory for output files
        normalize_intensities: Whether to normalize intensities before comparison
        min_peaks: Minimum number of peaks required for comparison
        ionization_filter: Filter by ionization mode (e.g., "EI", "ESI+")
        massspecgym_splits: Which MassSpecGym splits to use (train/val/test)
    """
    
    similarity_threshold: float = 0.5
    mz_tolerance: float = 0.5
    metrics: List[str] = field(default_factory=lambda: [
        "cosine", "dot_product", "spectral_entropy", "weighted_cosine"
    ])
    max_experimental_matches: int = 50
    output_dir: str = "./benchmark_results"
    normalize_intensities: bool = True
    min_peaks: int = 3
    ionization_filter: Optional[str] = None
    massspecgym_splits: List[str] = field(default_factory=lambda: ["train", "val", "test"])
    
    # Advanced options
    mz_power: float = 0.0  # For weighted cosine
    intensity_power: float = 0.5  # For weighted cosine
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        return asdict(self)
    
    def __post_init__(self):
        """Validate configuration."""
        if not 0 <= self.similarity_threshold <= 1:
            raise ValueError("similarity_threshold must be between 0 and 1")
        if self.mz_tolerance <= 0:
            raise ValueError("mz_tolerance must be positive")
        if self.max_experimental_matches < 1:
            raise ValueError("max_experimental_matches must be at least 1")


@dataclass
class ExperimentalMatch:
    """A single match between calculated and experimental spectrum."""
    
    experimental_id: str
    experimental_source: str
    experimental_smiles: Optional[str]
    experimental_formula: Optional[str]
    experimental_inchikey: Optional[str]
    experimental_adduct: Optional[str]
    experimental_peaks: int
    
    # Similarity metrics
    metrics: SimilarityMetrics
    
    # Quality assessment
    quality_tier: str = ""
    confidence: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "experimental_id": self.experimental_id,
            "experimental_source": self.experimental_source,
            "experimental_smiles": self.experimental_smiles,
            "experimental_formula": self.experimental_formula,
            "experimental_inchikey": self.experimental_inchikey,
            "experimental_adduct": self.experimental_adduct,
            "experimental_peaks": self.experimental_peaks,
            "metrics": self.metrics.to_dict(),
            "quality_tier": self.quality_tier,
            "confidence": self.confidence,
        }


@dataclass
class MoleculeBenchmark:
    """Benchmark results for a single molecule."""
    
    # Molecule identifiers
    smiles: str
    formula: str
    inchikey: str
    name: Optional[str] = None
    
    # Calculated spectrum info
    calculated_peaks: List[Tuple[float, float]] = field(default_factory=list)
    calculated_source: str = "QCxMS2"
    calculation_method: str = ""  # e.g., "EI", "CID"
    
    # Matching results
    experimental_matches: List[ExperimentalMatch] = field(default_factory=list)
    best_match: Optional[ExperimentalMatch] = None
    
    # Aggregate statistics
    best_cosine: float = 0.0
    best_spectral_entropy: float = 0.0
    mean_cosine: float = 0.0
    mean_spectral_entropy: float = 0.0
    num_experimental_found: int = 0
    num_above_threshold: int = 0
    
    # Metadata
    benchmark_time: float = 0.0  # Time to complete benchmark (seconds)
    search_method: str = ""  # How experimental spectra were found
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "smiles": self.smiles,
            "formula": self.formula,
            "inchikey": self.inchikey,
            "name": self.name,
            "calculated_peaks": self.calculated_peaks,
            "calculated_source": self.calculated_source,
            "calculation_method": self.calculation_method,
            "best_cosine": self.best_cosine,
            "best_spectral_entropy": self.best_spectral_entropy,
            "mean_cosine": self.mean_cosine,
            "mean_spectral_entropy": self.mean_spectral_entropy,
            "num_experimental_found": self.num_experimental_found,
            "num_above_threshold": self.num_above_threshold,
            "benchmark_time": self.benchmark_time,
            "search_method": self.search_method,
            "experimental_matches": [m.to_dict() for m in self.experimental_matches],
            "best_match": self.best_match.to_dict() if self.best_match else None,
        }


@dataclass
class BenchmarkResult:
    """
    Complete benchmark results including all molecules and summary statistics.
    """
    
    # Configuration
    config: BenchmarkConfig
    
    # Individual molecule results
    molecules: List[MoleculeBenchmark] = field(default_factory=list)
    
    # Aggregate statistics
    total_molecules: int = 0
    molecules_with_matches: int = 0
    molecules_above_threshold: int = 0
    
    # Metric statistics
    cosine_stats: Optional[StatisticalAnalysis] = None
    spectral_entropy_stats: Optional[StatisticalAnalysis] = None
    dot_product_stats: Optional[StatisticalAnalysis] = None
    
    # Metadata
    benchmark_name: str = ""
    benchmark_date: str = ""
    total_time: float = 0.0
    massspecgym_version: str = "NeurIPS 2024"
    qcxms2_version: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "benchmark_name": self.benchmark_name,
            "benchmark_date": self.benchmark_date,
            "total_time": self.total_time,
            "config": self.config.to_dict(),
            "summary": {
                "total_molecules": self.total_molecules,
                "molecules_with_matches": self.molecules_with_matches,
                "molecules_above_threshold": self.molecules_above_threshold,
                "match_rate": self.molecules_with_matches / max(1, self.total_molecules),
                "success_rate": self.molecules_above_threshold / max(1, self.total_molecules),
            },
            "statistics": {
                "cosine": self.cosine_stats.to_dict() if self.cosine_stats else None,
                "spectral_entropy": self.spectral_entropy_stats.to_dict() if self.spectral_entropy_stats else None,
                "dot_product": self.dot_product_stats.to_dict() if self.dot_product_stats else None,
            },
            "molecules": [m.to_dict() for m in self.molecules],
            "massspecgym_version": self.massspecgym_version,
            "qcxms2_version": self.qcxms2_version,
        }


class QCxMS2Benchmark:
    """
    Main benchmark class for comparing QCxMS2 calculated spectra against MassSpecGym.
    
    Example usage:
        benchmark = QCxMS2Benchmark(config)
        benchmark.add_molecule(smiles="CCO", formula="C2H6O", ...)
        results = benchmark.run_benchmark()
        benchmark.generate_report()
        benchmark.export_results("results.json")
    """
    
    def __init__(self, config: Optional[BenchmarkConfig] = None):
        """
        Initialize benchmark suite.
        
        Args:
            config: Benchmark configuration (uses defaults if not provided)
        """
        self.config = config or BenchmarkConfig()
        self._molecules: List[Dict[str, Any]] = []
        self._massspecgym_connector = None
        self._metrics_calculator = MetricsCalculator(
            mz_tolerance=self.config.mz_tolerance,
            mz_power=self.config.mz_power,
            intensity_power=self.config.intensity_power,
        )
        self._results: Optional[BenchmarkResult] = None
        
        # Ensure output directory exists
        Path(self.config.output_dir).mkdir(parents=True, exist_ok=True)
        
        logger.info(f"QCxMS2 Benchmark initialized with config: {self.config}")
    
    def _ensure_massspecgym_loaded(self) -> bool:
        """Ensure MassSpecGym connector is loaded."""
        if self._massspecgym_connector is not None:
            return True
        
        try:
            # Try multiple import strategies
            MassSpecGymConnector = None
            DatabaseConfig = None
            
            # Strategy 1: Direct relative import (when used as part of database_integration package)
            try:
                from ..connectors.massspecgym import MassSpecGymConnector
                from ..config import DatabaseConfig
            except (ImportError, ValueError):
                pass
            
            # Strategy 2: Absolute import (when database_integration is in sys.path)
            if MassSpecGymConnector is None:
                try:
                    from database_integration.connectors.massspecgym import MassSpecGymConnector
                    from database_integration.config import DatabaseConfig
                except ImportError:
                    pass
            
            # Strategy 3: Add parent to path and import directly
            if MassSpecGymConnector is None:
                import sys
                from pathlib import Path
                parent_dir = str(Path(__file__).parent.parent)
                if parent_dir not in sys.path:
                    sys.path.insert(0, parent_dir)
                from connectors.massspecgym import MassSpecGymConnector
                from config import DatabaseConfig
            
            if MassSpecGymConnector is None or DatabaseConfig is None:
                raise ImportError("Could not import MassSpecGymConnector")
            
            # Create connector with specified splits
            config = DatabaseConfig.for_massspecgym()
            config.additional_params["splits"] = self.config.massspecgym_splits
            
            self._massspecgym_connector = MassSpecGymConnector(config)
            
            # Test connection (this loads the dataset)
            if not self._massspecgym_connector.test_connection():
                logger.error("Failed to connect to MassSpecGym")
                return False
            
            logger.info("MassSpecGym connector loaded successfully")
            return True
            
        except ImportError as e:
            logger.error(f"Failed to import MassSpecGym connector: {e}")
            return False
        except Exception as e:
            logger.error(f"Failed to load MassSpecGym: {e}")
            return False
    
    def add_molecule(
        self,
        smiles: str,
        formula: str,
        inchikey: str,
        calculated_spectrum: Union[Dict[str, Any], List[Tuple[float, float]]],
        name: Optional[str] = None,
        calculation_method: str = "EI",
    ) -> None:
        """
        Add a molecule to the benchmark queue.
        
        Args:
            smiles: SMILES string
            formula: Molecular formula (e.g., "C6H12O6")
            inchikey: InChI key (full 27 chars or 14-char 2D key)
            calculated_spectrum: Either:
                - Dict with "peaks" key: {"peaks": [(mz, intensity), ...]}
                - List of (mz, intensity) tuples
            name: Optional molecule name
            calculation_method: QCxMS2 method used ("EI", "CID", etc.)
        """
        # Normalize spectrum format
        if isinstance(calculated_spectrum, dict):
            peaks = calculated_spectrum.get("peaks", [])
        else:
            peaks = calculated_spectrum
        
        # Validate peaks
        if not peaks:
            logger.warning(f"No peaks provided for {smiles}, skipping")
            return
        
        if len(peaks) < self.config.min_peaks:
            logger.warning(
                f"Too few peaks ({len(peaks)}) for {smiles}, "
                f"minimum is {self.config.min_peaks}"
            )
            return
        
        molecule = {
            "smiles": smiles,
            "formula": formula,
            "inchikey": inchikey,
            "calculated_peaks": peaks,
            "name": name,
            "calculation_method": calculation_method,
        }
        
        self._molecules.append(molecule)
        logger.debug(f"Added molecule: {smiles} ({formula}) with {len(peaks)} peaks")
    
    def add_molecules_from_file(
        self,
        filepath: Union[str, Path],
        spectrum_dir: Optional[Union[str, Path]] = None,
    ) -> int:
        """
        Add molecules from a CSV or JSON file.
        
        CSV format: smiles,formula,inchikey,name,spectrum_file
        JSON format: [{"smiles": ..., "formula": ..., "peaks": [...]}]
        
        Args:
            filepath: Path to input file
            spectrum_dir: Directory containing spectrum files (for CSV)
        
        Returns:
            Number of molecules added
        """
        filepath = Path(filepath)
        
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        added = 0
        
        if filepath.suffix.lower() == ".json":
            with open(filepath, "r") as f:
                data = json.load(f)
            
            for mol in data:
                try:
                    self.add_molecule(
                        smiles=mol["smiles"],
                        formula=mol["formula"],
                        inchikey=mol.get("inchikey", ""),
                        calculated_spectrum=mol.get("peaks", mol.get("calculated_spectrum", [])),
                        name=mol.get("name"),
                        calculation_method=mol.get("calculation_method", "EI"),
                    )
                    added += 1
                except (KeyError, ValueError) as e:
                    logger.warning(f"Failed to add molecule from JSON: {e}")
        
        elif filepath.suffix.lower() == ".csv":
            with open(filepath, "r", newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        # Load spectrum from file if specified
                        peaks = []
                        spectrum_file = row.get("spectrum_file", "")
                        if spectrum_file and spectrum_dir:
                            peaks = self._load_spectrum_file(
                                Path(spectrum_dir) / spectrum_file
                            )
                        
                        self.add_molecule(
                            smiles=row["smiles"],
                            formula=row["formula"],
                            inchikey=row.get("inchikey", ""),
                            calculated_spectrum=peaks,
                            name=row.get("name"),
                            calculation_method=row.get("calculation_method", "EI"),
                        )
                        added += 1
                    except (KeyError, ValueError) as e:
                        logger.warning(f"Failed to add molecule from CSV: {e}")
        
        else:
            raise ValueError(f"Unsupported file format: {filepath.suffix}")
        
        logger.info(f"Added {added} molecules from {filepath}")
        return added
    
    def _load_spectrum_file(self, filepath: Path) -> List[Tuple[float, float]]:
        """Load spectrum from file (simple two-column format)."""
        peaks = []
        if not filepath.exists():
            return peaks
        
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        peaks.append((mz, intensity))
                    except ValueError:
                        continue
        
        return peaks
    
    def run_benchmark(
        self,
        progress_callback: Optional[Callable[[int, int, str], None]] = None,
    ) -> BenchmarkResult:
        """
        Run benchmark on all added molecules.
        
        Compares each calculated spectrum against experimental spectra
        from MassSpecGym, computing similarity metrics.
        
        Args:
            progress_callback: Optional callback(current, total, smiles) for progress
        
        Returns:
            BenchmarkResult with all results and statistics
        """
        if not self._molecules:
            raise ValueError("No molecules added to benchmark. Use add_molecule() first.")
        
        if not self._ensure_massspecgym_loaded():
            raise RuntimeError(
                "Failed to load MassSpecGym. Ensure 'datasets' package is installed."
            )
        
        start_time = time.time()
        results: List[MoleculeBenchmark] = []
        
        total = len(self._molecules)
        logger.info(f"Starting benchmark on {total} molecules...")
        
        for i, mol in enumerate(self._molecules):
            if progress_callback:
                progress_callback(i + 1, total, mol["smiles"])
            
            logger.debug(f"Benchmarking {i+1}/{total}: {mol['smiles']}")
            
            try:
                result = self._benchmark_molecule(mol)
                results.append(result)
            except Exception as e:
                logger.error(f"Failed to benchmark {mol['smiles']}: {e}")
                # Add failed result
                results.append(MoleculeBenchmark(
                    smiles=mol["smiles"],
                    formula=mol["formula"],
                    inchikey=mol["inchikey"],
                    name=mol.get("name"),
                    calculated_peaks=mol["calculated_peaks"],
                ))
        
        total_time = time.time() - start_time
        
        # Compute aggregate statistics
        benchmark_result = self._compute_aggregate_statistics(results, total_time)
        
        self._results = benchmark_result
        
        logger.info(
            f"Benchmark complete: {benchmark_result.molecules_with_matches}/{total} "
            f"molecules matched, {benchmark_result.molecules_above_threshold} above threshold"
        )
        
        return benchmark_result
    
    def _benchmark_molecule(self, mol: Dict[str, Any]) -> MoleculeBenchmark:
        """Benchmark a single molecule against MassSpecGym."""
        mol_start = time.time()
        
        # Convert calculated peaks to Peak objects
        calc_peaks = [Peak(mz=mz, intensity=intensity) for mz, intensity in mol["calculated_peaks"]]
        
        # Normalize intensities if configured
        if self.config.normalize_intensities:
            max_intensity = max(p.intensity for p in calc_peaks) if calc_peaks else 1.0
            if max_intensity > 0:
                calc_peaks = [Peak(mz=p.mz, intensity=p.intensity / max_intensity * 100) for p in calc_peaks]
        
        # Search for experimental spectra
        experimental_spectra = self._find_experimental_spectra(mol)
        
        # Compare against each experimental spectrum
        matches: List[ExperimentalMatch] = []
        
        for exp_spectrum in experimental_spectra[:self.config.max_experimental_matches]:
            try:
                # Convert experimental peaks
                exp_peaks = [Peak(mz=p.mz, intensity=p.intensity) for p in exp_spectrum.peaks]
                
                # Normalize experimental intensities
                if self.config.normalize_intensities and exp_peaks:
                    max_exp = max(p.intensity for p in exp_peaks)
                    if max_exp > 0:
                        exp_peaks = [Peak(mz=p.mz, intensity=p.intensity / max_exp * 100) for p in exp_peaks]
                
                # Calculate metrics
                metrics = self._metrics_calculator.calculate_all_metrics(calc_peaks, exp_peaks)
                
                # Assess quality
                quality = assess_match_quality(metrics)
                
                match = ExperimentalMatch(
                    experimental_id=exp_spectrum.spectrum_id,
                    experimental_source=exp_spectrum.source,
                    experimental_smiles=exp_spectrum.smiles,
                    experimental_formula=exp_spectrum.formula,
                    experimental_inchikey=exp_spectrum.inchi_key,
                    experimental_adduct=exp_spectrum.precursor_type,
                    experimental_peaks=len(exp_spectrum.peaks),
                    metrics=metrics,
                    quality_tier=quality["quality_tier"],
                    confidence=quality["confidence"],
                )
                
                matches.append(match)
                
            except Exception as e:
                logger.debug(f"Failed to compare with {exp_spectrum.spectrum_id}: {e}")
                continue
        
        # Sort matches by cosine similarity
        matches.sort(key=lambda m: m.metrics.cosine, reverse=True)
        
        # Compute statistics for this molecule
        best_match = matches[0] if matches else None
        
        cosine_scores = [m.metrics.cosine for m in matches]
        entropy_scores = [m.metrics.spectral_entropy for m in matches]
        
        result = MoleculeBenchmark(
            smiles=mol["smiles"],
            formula=mol["formula"],
            inchikey=mol["inchikey"],
            name=mol.get("name"),
            calculated_peaks=mol["calculated_peaks"],
            calculation_method=mol.get("calculation_method", ""),
            experimental_matches=matches,
            best_match=best_match,
            best_cosine=max(cosine_scores) if cosine_scores else 0.0,
            best_spectral_entropy=max(entropy_scores) if entropy_scores else 0.0,
            mean_cosine=sum(cosine_scores) / len(cosine_scores) if cosine_scores else 0.0,
            mean_spectral_entropy=sum(entropy_scores) / len(entropy_scores) if entropy_scores else 0.0,
            num_experimental_found=len(matches),
            num_above_threshold=sum(1 for m in matches if m.metrics.cosine >= self.config.similarity_threshold),
            benchmark_time=time.time() - mol_start,
            search_method="inchikey" if mol["inchikey"] else "formula",
        )
        
        return result
    
    def _find_experimental_spectra(self, mol: Dict[str, Any]) -> List:
        """Find experimental spectra for a molecule from MassSpecGym."""
        # Import with multiple strategies for robustness
        try:
            from ..search_criteria import SpectrumSearchCriteria
        except (ImportError, ValueError):
            try:
                from database_integration.search_criteria import SpectrumSearchCriteria
            except ImportError:
                from search_criteria import SpectrumSearchCriteria
        
        # Try InChI key first (most specific)
        if mol["inchikey"]:
            criteria = SpectrumSearchCriteria(
                inchi_key=mol["inchikey"],
                max_results=self.config.max_experimental_matches,
            )
            result = self._massspecgym_connector.search_spectrum(criteria)
            if result.success and result.spectra:
                return result.spectra
        
        # Fall back to formula search
        if mol["formula"]:
            criteria = SpectrumSearchCriteria(
                formula=mol["formula"],
                max_results=self.config.max_experimental_matches * 2,  # Get more for filtering
            )
            result = self._massspecgym_connector.search_spectrum(criteria)
            if result.success and result.spectra:
                return result.spectra
        
        # Fall back to SMILES search
        if mol["smiles"]:
            criteria = SpectrumSearchCriteria(
                smiles=mol["smiles"],
                max_results=self.config.max_experimental_matches,
            )
            result = self._massspecgym_connector.search_spectrum(criteria)
            if result.success and result.spectra:
                return result.spectra
        
        return []
    
    def _compute_aggregate_statistics(
        self,
        results: List[MoleculeBenchmark],
        total_time: float,
    ) -> BenchmarkResult:
        """Compute aggregate statistics from all molecule benchmarks."""
        # Collect best scores from each molecule
        best_cosines = [r.best_cosine for r in results if r.best_cosine > 0]
        best_entropies = [r.best_spectral_entropy for r in results if r.best_spectral_entropy > 0]
        
        # Collect all match cosines for dot product stats
        all_cosines = []
        all_entropies = []
        all_dot_products = []
        
        for r in results:
            for m in r.experimental_matches:
                all_cosines.append(m.metrics.cosine)
                all_entropies.append(m.metrics.spectral_entropy)
                all_dot_products.append(m.metrics.dot_product)
        
        return BenchmarkResult(
            config=self.config,
            molecules=results,
            total_molecules=len(results),
            molecules_with_matches=sum(1 for r in results if r.num_experimental_found > 0),
            molecules_above_threshold=sum(1 for r in results if r.best_cosine >= self.config.similarity_threshold),
            cosine_stats=compute_statistical_analysis(best_cosines) if best_cosines else None,
            spectral_entropy_stats=compute_statistical_analysis(best_entropies) if best_entropies else None,
            dot_product_stats=compute_statistical_analysis(all_dot_products) if all_dot_products else None,
            benchmark_name=f"QCxMS2_benchmark_{datetime.now().strftime('%Y%m%d')}",
            benchmark_date=datetime.now().isoformat(),
            total_time=total_time,
        )
    
    def generate_report(self, output_path: Optional[Union[str, Path]] = None) -> str:
        """
        Generate a human-readable benchmark report in markdown format.
        
        Args:
            output_path: Optional path to save report (defaults to output_dir/report.md)
        
        Returns:
            Report as markdown string
        """
        if self._results is None:
            raise ValueError("No benchmark results. Run run_benchmark() first.")
        
        results = self._results
        
        lines = [
            "# QCxMS2 Benchmark Report",
            "",
            f"**Date:** {results.benchmark_date}",
            f"**Total Time:** {results.total_time:.2f} seconds",
            "",
            "## Summary",
            "",
            f"- **Total Molecules:** {results.total_molecules}",
            f"- **Molecules with Matches:** {results.molecules_with_matches} ({results.molecules_with_matches/max(1, results.total_molecules)*100:.1f}%)",
            f"- **Molecules Above Threshold ({self.config.similarity_threshold}):** {results.molecules_above_threshold} ({results.molecules_above_threshold/max(1, results.total_molecules)*100:.1f}%)",
            "",
            "## Configuration",
            "",
            f"- Similarity Threshold: {self.config.similarity_threshold}",
            f"- m/z Tolerance: {self.config.mz_tolerance} Da",
            f"- Metrics: {', '.join(self.config.metrics)}",
            f"- MassSpecGym Splits: {', '.join(self.config.massspecgym_splits)}",
            "",
        ]
        
        # Cosine similarity statistics
        if results.cosine_stats:
            stats = results.cosine_stats
            lines.extend([
                "## Cosine Similarity Statistics (Best Match per Molecule)",
                "",
                "| Statistic | Value |",
                "|-----------|-------|",
                f"| Count | {stats.count} |",
                f"| Mean | {stats.mean:.4f} |",
                f"| Std Dev | {stats.std:.4f} |",
                f"| Median | {stats.median:.4f} |",
                f"| Min | {stats.min_value:.4f} |",
                f"| Max | {stats.max_value:.4f} |",
                f"| 25th Percentile | {stats.percentile_25:.4f} |",
                f"| 75th Percentile | {stats.percentile_75:.4f} |",
                f"| 90th Percentile | {stats.percentile_90:.4f} |",
                f"| 95% CI | [{stats.ci_lower:.4f}, {stats.ci_upper:.4f}] |",
                "",
                "### Quality Distribution",
                "",
                f"- Above 0.9: {stats.above_threshold_90:.1f}%",
                f"- Above 0.7: {stats.above_threshold_70:.1f}%",
                f"- Above 0.5: {stats.above_threshold_50:.1f}%",
                "",
            ])
        
        # Top matches table
        lines.extend([
            "## Top 20 Best Matches",
            "",
            "| Molecule | Formula | Best Cosine | Best Entropy | Matches Found |",
            "|----------|---------|-------------|--------------|---------------|",
        ])
        
        # Sort by best cosine
        sorted_mols = sorted(results.molecules, key=lambda m: m.best_cosine, reverse=True)
        
        for mol in sorted_mols[:20]:
            name = mol.name or mol.smiles[:30]
            lines.append(
                f"| {name} | {mol.formula} | {mol.best_cosine:.4f} | "
                f"{mol.best_spectral_entropy:.4f} | {mol.num_experimental_found} |"
            )
        
        lines.extend(["", ""])
        
        # Poor matches (for investigation)
        poor_matches = [m for m in results.molecules if 0 < m.best_cosine < 0.5]
        if poor_matches:
            lines.extend([
                "## Molecules with Poor Matches (cosine < 0.5)",
                "",
                "| Molecule | Formula | Best Cosine | Matches Found |",
                "|----------|---------|-------------|---------------|",
            ])
            
            for mol in sorted(poor_matches, key=lambda m: m.best_cosine)[:10]:
                name = mol.name or mol.smiles[:30]
                lines.append(
                    f"| {name} | {mol.formula} | {mol.best_cosine:.4f} | {mol.num_experimental_found} |"
                )
            
            lines.extend(["", ""])
        
        # No matches
        no_matches = [m for m in results.molecules if m.num_experimental_found == 0]
        if no_matches:
            lines.extend([
                "## Molecules with No Experimental Matches",
                "",
            ])
            for mol in no_matches[:20]:
                lines.append(f"- {mol.smiles} ({mol.formula})")
            
            if len(no_matches) > 20:
                lines.append(f"- ... and {len(no_matches) - 20} more")
            
            lines.extend(["", ""])
        
        # Citation
        lines.extend([
            "## Data Sources",
            "",
            "- **MassSpecGym:** Bushuiev et al. MassSpecGym: A benchmark for the discovery and identification of molecules. NeurIPS 2024 (Spotlight).",
            "- **QCxMS2:** Grimme et al. Quantum Chemical Mass Spectrometry calculations.",
            "",
        ])
        
        report = "\n".join(lines)
        
        # Save report
        if output_path is None:
            output_path = Path(self.config.output_dir) / "benchmark_report.md"
        else:
            output_path = Path(output_path)
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(report)
        
        logger.info(f"Report saved to {output_path}")
        
        return report
    
    def export_results(
        self,
        filepath: Union[str, Path],
        format: str = "auto",
    ) -> None:
        """
        Export benchmark results to file.
        
        Args:
            filepath: Output file path
            format: Output format ("json", "csv", or "auto" to detect from extension)
        """
        if self._results is None:
            raise ValueError("No benchmark results. Run run_benchmark() first.")
        
        filepath = Path(filepath)
        
        # Auto-detect format
        if format == "auto":
            format = filepath.suffix.lower().lstrip(".")
        
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        if format == "json":
            self._export_json(filepath)
        elif format == "csv":
            self._export_csv(filepath)
        else:
            raise ValueError(f"Unsupported export format: {format}")
        
        logger.info(f"Results exported to {filepath}")
    
    def _export_json(self, filepath: Path) -> None:
        """Export results as JSON."""
        with open(filepath, "w") as f:
            json.dump(self._results.to_dict(), f, indent=2)
    
    def _export_csv(self, filepath: Path) -> None:
        """Export results as CSV (one row per molecule)."""
        fieldnames = [
            "smiles", "formula", "inchikey", "name",
            "best_cosine", "best_spectral_entropy", "mean_cosine",
            "num_experimental_found", "num_above_threshold",
            "best_match_id", "best_match_source",
            "calculated_peaks_count", "calculation_method",
        ]
        
        with open(filepath, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for mol in self._results.molecules:
                row = {
                    "smiles": mol.smiles,
                    "formula": mol.formula,
                    "inchikey": mol.inchikey,
                    "name": mol.name or "",
                    "best_cosine": mol.best_cosine,
                    "best_spectral_entropy": mol.best_spectral_entropy,
                    "mean_cosine": mol.mean_cosine,
                    "num_experimental_found": mol.num_experimental_found,
                    "num_above_threshold": mol.num_above_threshold,
                    "best_match_id": mol.best_match.experimental_id if mol.best_match else "",
                    "best_match_source": mol.best_match.experimental_source if mol.best_match else "",
                    "calculated_peaks_count": len(mol.calculated_peaks),
                    "calculation_method": mol.calculation_method,
                }
                writer.writerow(row)
    
    def get_results(self) -> Optional[BenchmarkResult]:
        """Get benchmark results (None if benchmark hasn't been run)."""
        return self._results
    
    def clear_molecules(self) -> None:
        """Clear all added molecules (to start fresh)."""
        self._molecules = []
        self._results = None
        logger.info("Cleared all molecules from benchmark")
    
    @property
    def num_molecules(self) -> int:
        """Number of molecules added to benchmark."""
        return len(self._molecules)
