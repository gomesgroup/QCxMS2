"""
MassBank connector for accessing MassBank EU database.

MassBank is an open-access mass spectra database with 40K+ spectra.
API Documentation: https://github.com/MassBank/MassBank-web/blob/main/Documentation/MassBankRecordFormat.md
"""

import logging
from typing import List, Optional

from ..base import MSDatabaseConnector
from ..config import DatabaseConfig
from ..data_classes import ExperimentalSpectrum, MSSearchResult, Peak
from ..enums import IonizationMode, SpectrumType
from ..search_criteria import SpectrumSearchCriteria

logger = logging.getLogger(__name__)


class MassBankConnector(MSDatabaseConnector):
    """
    Connector for MassBank EU database.

    Features:
    - REST API access to 40K+ mass spectra
    - Search by formula, exact mass, InChI, peaks
    - Supports EI, ESI+, ESI-, and other ionization modes
    - Free access, no API key required
    """

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize MassBank connector.

        Args:
            config: Optional custom configuration
        """
        if config is None:
            config = DatabaseConfig.for_massbank()
        super().__init__(config)

    def test_connection(self) -> bool:
        """Test connection to MassBank API."""
        try:
            # Try to get database info endpoint
            result = self._make_request(f"{self.config.base_url}/info")
            return result.get("success", False)
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

        # Determine search endpoint based on criteria
        if criteria.formula:
            return self._search_by_formula(criteria)
        elif criteria.exact_mass is not None:
            return self._search_by_exact_mass(criteria)
        elif criteria.inchi:
            return self._search_by_inchi(criteria)
        elif criteria.inchi_key:
            return self._search_by_inchi_key(criteria)
        elif criteria.peaks:
            return self._search_by_peaks(criteria)
        else:
            return MSSearchResult(
                success=False,
                error_message="MassBank requires formula, exact mass, InChI, InChI Key, or peaks for search",
                source="MassBank EU",
            )

    def _search_by_formula(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """Search by molecular formula."""
        url = f"{self.config.base_url}/searchSpectrum"

        params = {
            "formula": criteria.formula,
            "tolerance": criteria.exact_mass_tolerance,
        }

        if criteria.ionization_mode:
            params["ionMode"] = self._convert_ionization_mode(criteria.ionization_mode)

        if criteria.max_results:
            params["maxResults"] = criteria.max_results

        response = self._make_request(url, params=params)

        if not response["success"]:
            return MSSearchResult(
                success=False,
                error_message=response.get("error_message", "Unknown error"),
                source="MassBank EU",
            )

        # Parse response and extract spectra
        spectra = self._parse_search_results(response["data"])

        return MSSearchResult(
            success=True,
            spectra=spectra,
            total_results=len(spectra),
            source="MassBank EU",
            citation=self.get_citation(),
            request_time=response.get("request_time", 0.0),
            cached=response.get("cached", False),
        )

    def _search_by_exact_mass(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """Search by exact mass."""
        url = f"{self.config.base_url}/searchSpectrum"

        params = {
            "exactMass": criteria.exact_mass,
            "tolerance": criteria.exact_mass_tolerance,
        }

        if criteria.ionization_mode:
            params["ionMode"] = self._convert_ionization_mode(criteria.ionization_mode)

        if criteria.max_results:
            params["maxResults"] = criteria.max_results

        response = self._make_request(url, params=params)

        if not response["success"]:
            return MSSearchResult(
                success=False,
                error_message=response.get("error_message", "Unknown error"),
                source="MassBank EU",
            )

        spectra = self._parse_search_results(response["data"])

        return MSSearchResult(
            success=True,
            spectra=spectra,
            total_results=len(spectra),
            source="MassBank EU",
            citation=self.get_citation(),
            request_time=response.get("request_time", 0.0),
            cached=response.get("cached", False),
        )

    def _search_by_inchi(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """Search by InChI."""
        url = f"{self.config.base_url}/searchSpectrum"

        params = {"inchi": criteria.inchi}

        if criteria.max_results:
            params["maxResults"] = criteria.max_results

        response = self._make_request(url, params=params)

        if not response["success"]:
            return MSSearchResult(
                success=False,
                error_message=response.get("error_message", "Unknown error"),
                source="MassBank EU",
            )

        spectra = self._parse_search_results(response["data"])

        return MSSearchResult(
            success=True,
            spectra=spectra,
            total_results=len(spectra),
            source="MassBank EU",
            citation=self.get_citation(),
            request_time=response.get("request_time", 0.0),
            cached=response.get("cached", False),
        )

    def _search_by_inchi_key(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """Search by InChI Key."""
        url = f"{self.config.base_url}/searchSpectrum"

        params = {"inchikey": criteria.inchi_key}

        if criteria.max_results:
            params["maxResults"] = criteria.max_results

        response = self._make_request(url, params=params)

        if not response["success"]:
            return MSSearchResult(
                success=False,
                error_message=response.get("error_message", "Unknown error"),
                source="MassBank EU",
            )

        spectra = self._parse_search_results(response["data"])

        return MSSearchResult(
            success=True,
            spectra=spectra,
            total_results=len(spectra),
            source="MassBank EU",
            citation=self.get_citation(),
            request_time=response.get("request_time", 0.0),
            cached=response.get("cached", False),
        )

    def _search_by_peaks(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """Search by peak list (similarity search)."""
        url = f"{self.config.base_url}/searchSpectrum"

        # Format peaks as "mz1:intensity1 mz2:intensity2 ..."
        peak_str = " ".join([f"{mz}:{intensity}" for mz, intensity in criteria.peaks])

        params = {
            "peaks": peak_str,
            "tolerance": criteria.mz_tolerance,
            "cutoff": criteria.similarity_threshold * 100,  # Convert to percentage
        }

        if criteria.ionization_mode:
            params["ionMode"] = self._convert_ionization_mode(criteria.ionization_mode)

        if criteria.max_results:
            params["maxResults"] = criteria.max_results

        response = self._make_request(url, params=params)

        if not response["success"]:
            return MSSearchResult(
                success=False,
                error_message=response.get("error_message", "Unknown error"),
                source="MassBank EU",
            )

        spectra = self._parse_search_results(response["data"])

        return MSSearchResult(
            success=True,
            spectra=spectra,
            total_results=len(spectra),
            source="MassBank EU",
            citation=self.get_citation(),
            request_time=response.get("request_time", 0.0),
            cached=response.get("cached", False),
        )

    def get_spectrum_by_id(self, spectrum_id: str) -> Optional[ExperimentalSpectrum]:
        """
        Retrieve spectrum by MassBank accession ID.

        Args:
            spectrum_id: MassBank accession (e.g., "MSBNK-BAFG-CSL2328300")

        Returns:
            ExperimentalSpectrum or None if not found
        """
        url = f"{self.config.base_url}/getRecord"
        params = {"id": spectrum_id}

        response = self._make_request(url, params=params)

        if not response["success"]:
            logger.warning(f"Failed to retrieve spectrum {spectrum_id}: {response.get('error_message')}")
            return None

        return self._parse_massbank_record(response["data"])

    def _parse_search_results(self, data: any) -> List[ExperimentalSpectrum]:
        """Parse search results from MassBank API."""
        spectra = []

        if not isinstance(data, list):
            return spectra

        for record in data:
            try:
                spectrum = self._parse_massbank_record(record)
                if spectrum:
                    spectra.append(spectrum)
            except Exception as e:
                logger.warning(f"Failed to parse MassBank record: {e}")
                continue

        return spectra

    def _parse_massbank_record(self, record: dict) -> Optional[ExperimentalSpectrum]:
        """Parse a single MassBank record into ExperimentalSpectrum."""
        try:
            # Extract basic info
            spectrum_id = record.get("accession", "")
            name = record.get("compound", {}).get("names", [""])[0]
            formula = record.get("compound", {}).get("formula", "")
            exact_mass = record.get("compound", {}).get("mass", None)
            smiles = record.get("compound", {}).get("smiles", "")
            inchi = record.get("compound", {}).get("inchi", "")
            inchi_key = record.get("compound", {}).get("inchikey", "")

            # Extract acquisition info
            acquisition = record.get("acquisition", {})
            instrument = acquisition.get("instrument", "")
            instrument_type = acquisition.get("instrument_type", "")

            # Extract MS info
            ms_info = record.get("ms", {})
            ionization_mode_str = ms_info.get("ionization_mode", "")
            ionization_mode = self._parse_ionization_mode(ionization_mode_str)

            precursor_mz = ms_info.get("precursor_mz", None)
            precursor_type = ms_info.get("precursor_type", "")
            collision_energy = ms_info.get("collision_energy", "")

            # Extract peaks
            peaks_data = record.get("peaks", {}).get("peak_list", [])
            peaks = [Peak(mz=float(p[0]), intensity=float(p[1])) for p in peaks_data if len(p) >= 2]

            # Extract metadata
            authors = record.get("contributor", "").split(", ") if record.get("contributor") else []
            date_created = record.get("date", {}).get("created", "")

            url = f"https://massbank.eu/MassBank/RecordDisplay?id={spectrum_id}"

            spectrum = ExperimentalSpectrum(
                spectrum_id=spectrum_id,
                name=name,
                formula=formula,
                exact_mass=exact_mass,
                smiles=smiles,
                inchi=inchi,
                inchi_key=inchi_key,
                peaks=peaks,
                precursor_mz=precursor_mz,
                precursor_type=precursor_type,
                ionization_mode=ionization_mode,
                spectrum_type=SpectrumType.EXPERIMENTAL,
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
        """Convert IonizationMode enum to MassBank format."""
        mapping = {
            IonizationMode.EI: "EI",
            IonizationMode.ESI_POSITIVE: "ESI+",
            IonizationMode.ESI_NEGATIVE: "ESI-",
            IonizationMode.APCI_POSITIVE: "APCI+",
            IonizationMode.APCI_NEGATIVE: "APCI-",
            IonizationMode.CI: "CI",
            IonizationMode.CID: "CID",
        }
        return mapping.get(mode, mode.value)

    @staticmethod
    def _parse_ionization_mode(mode_str: str) -> Optional[IonizationMode]:
        """Parse ionization mode string to enum."""
        mode_str_upper = mode_str.upper().strip()

        mapping = {
            "EI": IonizationMode.EI,
            "ESI+": IonizationMode.ESI_POSITIVE,
            "ESI-": IonizationMode.ESI_NEGATIVE,
            "APCI+": IonizationMode.APCI_POSITIVE,
            "APCI-": IonizationMode.APCI_NEGATIVE,
            "CI": IonizationMode.CI,
            "CID": IonizationMode.CID,
            "MALDI": IonizationMode.MALDI,
        }

        return mapping.get(mode_str_upper, IonizationMode.OTHER)

    def get_citation(self) -> str:
        """Get citation for MassBank."""
        return (
            "Horai et al. MassBank: a public repository for sharing mass spectral data "
            "for life sciences. J. Mass Spectrom. 2010, 45(7):703-714. "
            "https://massbank.eu/"
        )
