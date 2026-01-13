"""
Data classes for mass spectrometry database integration.

Defines data structures for MS spectra, peaks, search results, and comparisons.
"""

import time
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from .enums import IonizationMode, SimilarityMetric, SpectrumType


@dataclass
class Peak:
    """A single peak in a mass spectrum."""

    mz: float  # Mass-to-charge ratio
    intensity: float  # Relative intensity (0-100 or 0-1000)
    annotation: Optional[str] = None  # Optional fragment annotation (e.g., "M+", "M-H2O")

    def __post_init__(self):
        """Validate peak data."""
        if self.mz < 0:
            raise ValueError(f"m/z must be non-negative, got {self.mz}")
        if self.intensity < 0:
            raise ValueError(f"Intensity must be non-negative, got {self.intensity}")

    def normalize_intensity(self, max_intensity: float = 100.0) -> "Peak":
        """Return a new Peak with normalized intensity."""
        if self.intensity == 0:
            return Peak(self.mz, 0.0, self.annotation)
        normalized = (self.intensity / max(self.intensity, 1e-10)) * max_intensity
        return Peak(self.mz, normalized, self.annotation)


@dataclass
class ExperimentalSpectrum:
    """Experimental mass spectrum from a database."""

    spectrum_id: str  # Database-specific identifier
    name: Optional[str] = None  # Compound name
    formula: Optional[str] = None  # Molecular formula
    exact_mass: Optional[float] = None  # Exact mass
    smiles: Optional[str] = None  # SMILES string
    inchi: Optional[str] = None  # InChI string
    inchi_key: Optional[str] = None  # InChI Key
    peaks: List[Peak] = field(default_factory=list)  # List of peaks
    precursor_mz: Optional[float] = None  # Precursor m/z (for MS2)
    precursor_type: Optional[str] = None  # e.g., "[M]+", "[M+H]+"
    ionization_mode: Optional[IonizationMode] = None  # Ionization method
    spectrum_type: SpectrumType = SpectrumType.EXPERIMENTAL
    collision_energy: Optional[str] = None  # Collision energy (e.g., "40 eV")
    instrument: Optional[str] = None  # Instrument type
    instrument_type: Optional[str] = None  # e.g., "GC-EI-TOF"
    source: str = ""  # Database source
    database_id: Optional[str] = None  # Original database ID
    url: Optional[str] = None  # URL to spectrum in database
    date_created: Optional[str] = None  # Creation date
    authors: List[str] = field(default_factory=list)  # Authors
    reference: Optional[str] = None  # Citation/reference
    doi: Optional[str] = None  # DOI
    comments: Optional[str] = None  # Additional comments
    metadata: Dict[str, Any] = field(default_factory=dict)  # Extra metadata

    def __post_init__(self):
        """Validate and process spectrum data."""
        # Sort peaks by m/z
        self.peaks = sorted(self.peaks, key=lambda p: p.mz)

        # Calculate exact mass from formula if not provided
        if self.exact_mass is None and self.formula:
            try:
                self.exact_mass = self._calculate_exact_mass(self.formula)
            except Exception:
                pass  # Silently ignore formula parsing errors

    @staticmethod
    def _calculate_exact_mass(formula: str) -> float:
        """Calculate exact mass from molecular formula (basic implementation)."""
        # Atomic masses for common elements
        atomic_masses = {
            "H": 1.007825,
            "C": 12.0000,
            "N": 14.003074,
            "O": 15.994915,
            "F": 18.998403,
            "P": 30.973762,
            "S": 31.972071,
            "Cl": 34.968853,
            "Br": 78.918338,
            "I": 126.904473,
        }

        import re

        mass = 0.0
        # Simple regex parser (element + optional count)
        for match in re.finditer(r"([A-Z][a-z]?)(\d*)", formula):
            element = match.group(1)
            count = int(match.group(2)) if match.group(2) else 1
            if element in atomic_masses:
                mass += atomic_masses[element] * count

        return mass

    def get_max_intensity(self) -> float:
        """Get the maximum intensity in the spectrum."""
        if not self.peaks:
            return 0.0
        return max(p.intensity for p in self.peaks)

    def normalize(self, max_intensity: float = 100.0) -> "ExperimentalSpectrum":
        """Return a new spectrum with normalized intensities."""
        if not self.peaks:
            return self

        current_max = self.get_max_intensity()
        if current_max == 0:
            return self

        normalized_peaks = [
            Peak(p.mz, (p.intensity / current_max) * max_intensity, p.annotation)
            for p in self.peaks
        ]

        # Create a copy with normalized peaks
        return ExperimentalSpectrum(
            spectrum_id=self.spectrum_id,
            name=self.name,
            formula=self.formula,
            exact_mass=self.exact_mass,
            smiles=self.smiles,
            inchi=self.inchi,
            inchi_key=self.inchi_key,
            peaks=normalized_peaks,
            precursor_mz=self.precursor_mz,
            precursor_type=self.precursor_type,
            ionization_mode=self.ionization_mode,
            spectrum_type=self.spectrum_type,
            collision_energy=self.collision_energy,
            instrument=self.instrument,
            instrument_type=self.instrument_type,
            source=self.source,
            database_id=self.database_id,
            url=self.url,
            date_created=self.date_created,
            authors=self.authors,
            reference=self.reference,
            doi=self.doi,
            comments=self.comments,
            metadata=self.metadata.copy(),
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert spectrum to dictionary."""
        result = asdict(self)
        if self.ionization_mode:
            result["ionization_mode"] = self.ionization_mode.value
        if self.spectrum_type:
            result["spectrum_type"] = self.spectrum_type.value
        return result


@dataclass
class SpectrumComparisonResult:
    """Result of comparing two mass spectra."""

    similarity_score: float  # Overall similarity (0-1)
    metric: SimilarityMetric  # Similarity metric used
    matched_peaks: int  # Number of matched peaks
    total_peaks_query: int  # Total peaks in query spectrum
    total_peaks_reference: int  # Total peaks in reference spectrum
    query_spectrum_id: str  # ID of query spectrum
    reference_spectrum_id: str  # ID of reference spectrum
    matched_peak_pairs: List[Tuple[Peak, Peak]] = field(
        default_factory=list
    )  # Pairs of matched peaks
    mz_tolerance: float = 0.5  # m/z tolerance used for matching
    metadata: Dict[str, Any] = field(default_factory=dict)  # Additional comparison data

    def __post_init__(self):
        """Validate comparison result."""
        if not 0 <= self.similarity_score <= 1:
            raise ValueError(f"Similarity score must be between 0 and 1, got {self.similarity_score}")


@dataclass
class MSSearchResult:
    """Result from a mass spectrum database search."""

    success: bool
    spectra: List[ExperimentalSpectrum] = field(default_factory=list)  # Matching spectra
    total_results: int = 0  # Total number of results
    error_message: Optional[str] = None
    request_time: float = field(default_factory=time.time)
    source: Optional[str] = None  # Database source
    citation: Optional[str] = None  # Citation for database
    metadata: Dict[str, Any] = field(default_factory=dict)
    cached: bool = False
    cache_path: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary."""
        return {
            "success": self.success,
            "spectra": [s.to_dict() for s in self.spectra],
            "total_results": self.total_results,
            "error_message": self.error_message,
            "request_time": self.request_time,
            "source": self.source,
            "citation": self.citation,
            "metadata": self.metadata,
            "cached": self.cached,
            "cache_path": self.cache_path,
        }


@dataclass
class BatchComparisonResult:
    """Result of comparing a query spectrum against multiple reference spectra."""

    query_spectrum_id: str
    comparisons: List[SpectrumComparisonResult] = field(default_factory=list)
    top_matches: List[SpectrumComparisonResult] = field(default_factory=list)  # Top N matches
    metric: SimilarityMetric = SimilarityMetric.COSINE
    mz_tolerance: float = 0.5

    def get_top_n(self, n: int = 10) -> List[SpectrumComparisonResult]:
        """Get top N matches sorted by similarity score."""
        sorted_comparisons = sorted(self.comparisons, key=lambda x: x.similarity_score, reverse=True)
        return sorted_comparisons[:n]

    def filter_by_threshold(self, threshold: float = 0.7) -> List[SpectrumComparisonResult]:
        """Filter comparisons by minimum similarity score."""
        return [c for c in self.comparisons if c.similarity_score >= threshold]
