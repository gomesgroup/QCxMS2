"""
Search criteria for mass spectrometry database queries.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .enums import IonizationMode, SearchType


@dataclass
class SpectrumSearchCriteria:
    """
    Criteria for searching mass spectrometry databases.

    Supports various search modes: identity, similarity, formula, exact mass, etc.
    """

    # Molecular identifiers
    name: Optional[str] = None  # Compound name
    formula: Optional[str] = None  # Molecular formula (e.g., "C6H12O6")
    smiles: Optional[str] = None  # SMILES string
    inchi: Optional[str] = None  # InChI string
    inchi_key: Optional[str] = None  # InChI Key
    cas: Optional[str] = None  # CAS registry number

    # Mass-related search
    exact_mass: Optional[float] = None  # Exact mass
    exact_mass_tolerance: float = 0.01  # Tolerance for exact mass (Da)
    precursor_mz: Optional[float] = None  # Precursor m/z
    precursor_mz_tolerance: float = 0.5  # Tolerance for precursor m/z (Da)

    # Spectrum-based search
    peaks: List[Tuple[float, float]] = field(default_factory=list)  # [(mz, intensity), ...]
    similarity_threshold: float = 0.7  # Minimum similarity score (0-1)
    mz_tolerance: float = 0.5  # m/z tolerance for peak matching (Da)

    # MS experimental conditions
    ionization_mode: Optional[IonizationMode] = None  # e.g., EI, ESI+, CID
    collision_energy: Optional[str] = None  # e.g., "40 eV"
    collision_energy_range: Optional[Tuple[float, float]] = None  # Min/max CE
    instrument_type: Optional[str] = None  # e.g., "GC-MS", "LC-MS/MS"

    # Search behavior
    search_type: SearchType = SearchType.FORMULA  # Type of search
    max_results: int = 100  # Maximum number of results to return
    offset: int = 0  # Offset for pagination

    # Filters
    mass_range: Optional[Tuple[float, float]] = None  # (min_mass, max_mass)
    include_predicted: bool = False  # Include theoretical/predicted spectra
    include_insilico: bool = False  # Include in silico spectra

    # Additional filters
    authors: Optional[List[str]] = None  # Filter by authors
    publication_year_range: Optional[Tuple[int, int]] = None  # (min_year, max_year)
    license: Optional[str] = None  # Filter by data license

    def validate(self) -> bool:
        """
        Validate search criteria.

        Returns:
            True if criteria are valid
        
        Raises:
            ValueError if criteria are invalid
        """
        # At least one search parameter must be provided
        has_identifier = any(
            [
                self.name,
                self.formula,
                self.smiles,
                self.inchi,
                self.inchi_key,
                self.cas,
                self.exact_mass is not None,
                self.precursor_mz is not None,
                bool(self.peaks),
            ]
        )

        if not has_identifier:
            raise ValueError(
                "At least one search parameter must be provided "
                "(name, formula, SMILES, InChI, mass, or peaks)"
            )

        # Validate tolerances
        if self.exact_mass_tolerance < 0:
            raise ValueError(f"exact_mass_tolerance must be non-negative, got {self.exact_mass_tolerance}")

        if self.precursor_mz_tolerance < 0:
            raise ValueError(
                f"precursor_mz_tolerance must be non-negative, got {self.precursor_mz_tolerance}"
            )

        if self.mz_tolerance < 0:
            raise ValueError(f"mz_tolerance must be non-negative, got {self.mz_tolerance}")

        # Validate similarity threshold
        if not 0 <= self.similarity_threshold <= 1:
            raise ValueError(
                f"similarity_threshold must be between 0 and 1, got {self.similarity_threshold}"
            )

        return True
