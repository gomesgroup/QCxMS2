"""
MoNA (MassBank of North America) connector.

MoNA is a free metabolomics mass spectra repository.
API Documentation: https://mona.fiehnlab.ucdavis.edu/rest/api-doc.html
"""

import logging
from typing import List, Optional

from ..base import MSDatabaseConnector
from ..config import DatabaseConfig
from ..data_classes import ExperimentalSpectrum, MSSearchResult, Peak
from ..enums import IonizationMode, SpectrumType
from ..search_criteria import SpectrumSearchCriteria

logger = logging.getLogger(__name__)


class MoNAConnector(MSDatabaseConnector):
    """
    Connector for MassBank of North America (MoNA).

    Features:
    - REST API access to metabolomics spectra
    - Search by formula, InChI, similarity
    - Supports various ionization modes
    - Free access, no API key required
    """

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """Initialize MoNA connector."""
        if config is None:
            config = DatabaseConfig.for_mona()
        super().__init__(config)

    def test_connection(self) -> bool:
        """Test connection to MoNA API."""
        try:
            # Try a simple search
            result = self._make_request(f"{self.config.base_url}/spectra", params={"size": 1})
            return result.get("success", False)
        except Exception as e:
            logger.error(f"MoNA connection test failed: {e}")
            return False

    def search_spectrum(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """
        Search MoNA for spectra (implementation similar to MassBank).
        
        Note: This is a simplified stub. Full implementation would parse
        MoNA's specific API format and response structure.
        """
        criteria.validate()

        # Build query string based on criteria
        query_parts = []
        if criteria.formula:
            query_parts.append(f'compound.molFormula=="{criteria.formula}"')
        if criteria.inchi_key:
            query_parts.append(f'compound.inchiKey=="{criteria.inchi_key}"')

        if not query_parts:
            return MSSearchResult(
                success=False,
                error_message="MoNA requires formula or InChI Key for search",
                source="MoNA",
            )

        query = " AND ".join(query_parts)
        url = f"{self.config.base_url}/spectra/search"

        params = {
            "query": query,
            "size": criteria.max_results,
        }

        response = self._make_request(url, params=params)

        if not response["success"]:
            return MSSearchResult(
                success=False,
                error_message=response.get("error_message", "Unknown error"),
                source="MoNA",
            )

        # Parse results (simplified)
        spectra = self._parse_search_results(response["data"])

        return MSSearchResult(
            success=True,
            spectra=spectra,
            total_results=len(spectra),
            source="MoNA",
            citation=self.get_citation(),
            request_time=response.get("request_time", 0.0),
        )

    def get_spectrum_by_id(self, spectrum_id: str) -> Optional[ExperimentalSpectrum]:
        """Retrieve spectrum by MoNA ID."""
        url = f"{self.config.base_url}/spectra/{spectrum_id}"
        response = self._make_request(url)

        if not response["success"]:
            return None

        return self._parse_mona_record(response["data"])

    def _parse_search_results(self, data: any) -> List[ExperimentalSpectrum]:
        """Parse MoNA search results."""
        spectra = []
        if isinstance(data, list):
            for record in data:
                try:
                    spectrum = self._parse_mona_record(record)
                    if spectrum:
                        spectra.append(spectrum)
                except Exception as e:
                    logger.warning(f"Failed to parse MoNA record: {e}")
        return spectra

    def _parse_mona_record(self, record: dict) -> Optional[ExperimentalSpectrum]:
        """Parse a MoNA record (simplified implementation)."""
        try:
            spectrum_id = record.get("id", "")
            compound = record.get("compound", [{}])[0]
            name = compound.get("names", [{}])[0].get("name", "")
            formula = compound.get("molFormula", "")
            inchi = compound.get("inchi", "")
            inchi_key = compound.get("inchiKey", "")

            # Parse peaks
            peaks_str = record.get("spectrum", "")
            peaks = []
            for line in peaks_str.split("\n"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    peaks.append(Peak(mz=float(parts[0]), intensity=float(parts[1])))

            return ExperimentalSpectrum(
                spectrum_id=spectrum_id,
                name=name,
                formula=formula,
                inchi=inchi,
                inchi_key=inchi_key,
                peaks=peaks,
                source="MoNA",
                database_id=spectrum_id,
                url=f"https://mona.fiehnlab.ucdavis.edu/spectra/display/{spectrum_id}",
            )
        except Exception as e:
            logger.error(f"Error parsing MoNA record: {e}")
            return None

    def get_citation(self) -> str:
        """Get citation for MoNA."""
        return "MassBank of North America (MoNA). https://mona.fiehnlab.ucdavis.edu/"
