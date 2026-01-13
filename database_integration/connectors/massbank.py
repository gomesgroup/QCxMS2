"""
MassBank connector for accessing MassBank EU database (MassBank 3 API).

MassBank is an open-access mass spectra database with 100K+ spectra.
API Documentation: https://massbank.eu/MassBank-api/ui/
"""

import logging
import time
from typing import Any, Dict, List, Optional

import requests

from ..base import MSDatabaseConnector
from ..config import DatabaseConfig
from ..data_classes import ExperimentalSpectrum, MSSearchResult, Peak
from ..enums import IonizationMode, SpectrumType
from ..search_criteria import SpectrumSearchCriteria

logger = logging.getLogger(__name__)

# MassBank 3 API base URL
MASSBANK3_API_URL = "https://massbank.eu/MassBank-api"


class MassBankConnector(MSDatabaseConnector):
    """
    Connector for MassBank EU database (MassBank 3 API).

    Features:
    - REST API access to 100K+ mass spectra
    - Search by formula, exact mass, InChI, InChI key, peaks
    - Supports EI, ESI+, ESI-, CI, and other ionization modes
    - Free access, no API key required
    
    Note: This connector uses the MassBank 3 API (2024+), not the legacy API.
    """

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize MassBank connector.

        Args:
            config: Optional custom configuration
        """
        if config is None:
            config = DatabaseConfig.for_massbank()
        
        # Override base_url to use MassBank 3 API
        config.base_url = MASSBANK3_API_URL
        super().__init__(config)

    def test_connection(self) -> bool:
        """Test connection to MassBank API."""
        try:
            # Try to get a simple record list (limit 1)
            response = requests.get(
                f"{self.config.base_url}/records",
                params={"formula": "C6H6"},  # Benzene as test
                timeout=self.config.timeout,
            )
            return response.status_code == 200 and isinstance(response.json(), list)
        except Exception as e:
            logger.error(f"MassBank connection test failed: {e}")
            return False

    def search_spectrum(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """
        Search MassBank for spectra matching criteria.

        Args:
            criteria: Search criteria

        Returns:
            MSSearchResult with matching spectra
        """
        criteria.validate()

        # Build query parameters
        params = {}
        
        if criteria.formula:
            params["formula"] = criteria.formula
        
        if criteria.inchi_key:
            params["inchi_key"] = criteria.inchi_key
        
        if criteria.inchi:
            params["inchi"] = criteria.inchi
        
        if criteria.exact_mass is not None:
            params["exact_mass"] = criteria.exact_mass
            params["mass_tolerance"] = criteria.exact_mass_tolerance
        
        if criteria.ionization_mode:
            params["ion_mode"] = self._convert_ionization_mode(criteria.ionization_mode)
        
        # Peak list search
        if criteria.peaks:
            # Format as "mz1,mz2,mz3" 
            peak_mzs = ",".join([str(mz) for mz, _ in criteria.peaks])
            params["peaks"] = peak_mzs
        
        if not params:
            return MSSearchResult(
                success=False,
                error_message="MassBank requires formula, exact mass, InChI, InChI Key, or peaks for search",
                source="MassBank EU",
            )

        # Make API request
        start_time = time.time()
        try:
            response = requests.get(
                f"{self.config.base_url}/records",
                params=params,
                timeout=self.config.timeout,
                headers={"User-Agent": self.config.user_agent},
            )
            response.raise_for_status()
            data = response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"MassBank API request failed: {e}")
            return MSSearchResult(
                success=False,
                error_message=str(e),
                source="MassBank EU",
            )
        except ValueError as e:
            logger.error(f"Failed to parse MassBank response: {e}")
            return MSSearchResult(
                success=False,
                error_message=f"Invalid JSON response: {e}",
                source="MassBank EU",
            )

        request_time = time.time() - start_time

        # Parse response
        if not isinstance(data, list):
            return MSSearchResult(
                success=False,
                error_message=f"Unexpected response type: {type(data)}",
                source="MassBank EU",
            )

        # Convert to ExperimentalSpectrum objects
        spectra = []
        for record in data:
            try:
                spectrum = self._parse_massbank_record(record)
                if spectrum:
                    spectra.append(spectrum)
            except Exception as e:
                logger.warning(f"Failed to parse MassBank record: {e}")
                continue

        # Apply max_results limit
        if criteria.max_results and len(spectra) > criteria.max_results:
            spectra = spectra[:criteria.max_results]

        return MSSearchResult(
            success=True,
            spectra=spectra,
            total_results=len(spectra),
            source="MassBank EU",
            citation=self.get_citation(),
            request_time=request_time,
            cached=False,
        )

    def get_spectrum_by_id(self, spectrum_id: str) -> Optional[ExperimentalSpectrum]:
        """
        Retrieve spectrum by MassBank accession ID.

        Args:
            spectrum_id: MassBank accession (e.g., "MSBNK-BAFG-CSL2328300")

        Returns:
            ExperimentalSpectrum or None if not found
        """
        try:
            response = requests.get(
                f"{self.config.base_url}/records/{spectrum_id}",
                timeout=self.config.timeout,
                headers={"User-Agent": self.config.user_agent},
            )
            response.raise_for_status()
            record = response.json()
            return self._parse_massbank_record(record)
        except Exception as e:
            logger.warning(f"Failed to retrieve spectrum {spectrum_id}: {e}")
            return None

    def _parse_massbank_record(self, record: Dict[str, Any]) -> Optional[ExperimentalSpectrum]:
        """Parse a MassBank 3 API record into ExperimentalSpectrum."""
        try:
            # Extract basic info
            spectrum_id = record.get("accession", "")
            
            compound = record.get("compound", {})
            names = compound.get("names", [])
            name = names[0] if names else ""
            formula = compound.get("formula", "")
            exact_mass = compound.get("mass", None)
            smiles = compound.get("smiles", "")
            inchi = compound.get("inchi", "")
            
            # Extract InChI key from links
            inchi_key = ""
            links = compound.get("link", [])
            for link in links:
                if link.get("database") == "INCHIKEY":
                    inchi_key = link.get("identifier", "")
                    break

            # Extract acquisition info
            acquisition = record.get("acquisition", {})
            instrument = acquisition.get("instrument", "")
            instrument_type = acquisition.get("instrument_type", "")
            
            ms_info = acquisition.get("mass_spectrometry", {})
            ion_mode_str = ms_info.get("ion_mode", "")
            ionization_mode = self._parse_ionization_mode(ion_mode_str, instrument_type)

            # Extract collision energy from subtags
            collision_energy = ""
            subtags = ms_info.get("subtags", [])
            for subtag in subtags:
                if subtag.get("subtag") in ["COLLISION_ENERGY", "IONIZATION_ENERGY"]:
                    collision_energy = subtag.get("value", "")
                    break

            # Extract peaks
            peak_data = record.get("peak", {}).get("peak", {})
            peak_values = peak_data.get("values", [])
            
            peaks = []
            for p in peak_values:
                mz = p.get("mz")
                intensity = p.get("intensity") or p.get("rel")  # Try both
                if mz is not None and intensity is not None:
                    peaks.append(Peak(mz=float(mz), intensity=float(intensity)))

            if not peaks:
                return None

            # Extract metadata
            authors = [a.get("name", "") for a in record.get("authors", [])]
            date_info = record.get("date", {})
            date_created = date_info.get("created", "")

            url = f"https://massbank.eu/MassBank/RecordDisplay?id={spectrum_id}"

            # Determine spectrum type from instrument
            spectrum_type = SpectrumType.EXPERIMENTAL
            if "MS2" in str(ms_info.get("ms_type", "")) or "MS/MS" in instrument_type:
                spectrum_type = SpectrumType.MS2

            spectrum = ExperimentalSpectrum(
                spectrum_id=spectrum_id,
                name=name,
                formula=formula,
                exact_mass=exact_mass,
                smiles=smiles,
                inchi=inchi,
                inchi_key=inchi_key,
                peaks=peaks,
                ionization_mode=ionization_mode,
                spectrum_type=spectrum_type,
                collision_energy=collision_energy,
                instrument=instrument,
                instrument_type=instrument_type,
                source="MassBank EU",
                database_id=spectrum_id,
                url=url,
                date_created=date_created,
                authors=authors,
            )

            return spectrum

        except Exception as e:
            logger.error(f"Error parsing MassBank record: {e}")
            return None

    @staticmethod
    def _convert_ionization_mode(mode: IonizationMode) -> str:
        """Convert IonizationMode enum to MassBank API format."""
        mapping = {
            IonizationMode.EI: "POSITIVE",  # EI is typically positive mode
            IonizationMode.ESI_POSITIVE: "POSITIVE",
            IonizationMode.ESI_NEGATIVE: "NEGATIVE",
            IonizationMode.APCI_POSITIVE: "POSITIVE",
            IonizationMode.APCI_NEGATIVE: "NEGATIVE",
            IonizationMode.CI: "POSITIVE",
            IonizationMode.CID: "POSITIVE",
        }
        return mapping.get(mode, "POSITIVE")

    @staticmethod
    def _parse_ionization_mode(ion_mode_str: str, instrument_type: str = "") -> Optional[IonizationMode]:
        """Parse ionization mode from MassBank fields."""
        ion_mode_upper = ion_mode_str.upper().strip()
        instr_upper = instrument_type.upper()

        # Check instrument type for EI
        if "EI" in instr_upper:
            return IonizationMode.EI
        elif "CI" in instr_upper and "APCI" not in instr_upper:
            return IonizationMode.CI
        elif "ESI" in instr_upper:
            if ion_mode_upper == "NEGATIVE":
                return IonizationMode.ESI_NEGATIVE
            return IonizationMode.ESI_POSITIVE
        elif "APCI" in instr_upper:
            if ion_mode_upper == "NEGATIVE":
                return IonizationMode.APCI_NEGATIVE
            return IonizationMode.APCI_POSITIVE

        # Fall back to ion mode string
        if ion_mode_upper == "POSITIVE":
            return IonizationMode.ESI_POSITIVE
        elif ion_mode_upper == "NEGATIVE":
            return IonizationMode.ESI_NEGATIVE

        return IonizationMode.OTHER

    def get_citation(self) -> str:
        """Get citation for MassBank."""
        return (
            "Horai et al. MassBank: a public repository for sharing mass spectral data "
            "for life sciences. J. Mass Spectrom. 2010, 45(7):703-714. "
            "https://massbank.eu/"
        )
