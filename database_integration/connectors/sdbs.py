"""
SDBS (Spectral Database for Organic Compounds) connector (stub).

Note: SDBS has a 50 queries/day limit and no official API.
This is a stub implementation requiring web scraping.
"""

import logging
from typing import Optional

from ..base import MSDatabaseConnector
from ..config import DatabaseConfig
from ..data_classes import ExperimentalSpectrum, MSSearchResult
from ..search_criteria import SpectrumSearchCriteria

logger = logging.getLogger(__name__)


class SDBSConnector(MSDatabaseConnector):
    """
    Connector for SDBS (Spectral Database for Organic Compounds).

    Note: SDBS has strict usage limits (50 queries/day) and no public API.
    This is a stub implementation. For production use, consider MassBank or MoNA.
    """

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """Initialize SDBS connector."""
        if config is None:
            config = DatabaseConfig.for_sdbs()
        super().__init__(config)
        logger.warning(
            "SDBS has a 50 queries/day limit and no official API. Consider using Mass Bank or MoNA."
        )

    def test_connection(self) -> bool:
        """Test connection to SDBS."""
        try:
            result = self._make_request("https://sdbs.db.aist.go.jp")
            return result.get("success", False)
        except Exception as e:
            logger.error(f"SDBS connection test failed: {e}")
            return False

    def search_spectrum(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """
        Search SDBS (stub implementation).

        Note: SDBS has a 50/day limit. This requires web scraping which is not fully implemented.
        """
        return MSSearchResult(
            success=False,
            error_message="SDBS connector not fully implemented (no public API, 50/day limit). Use MassBank or MoNA instead.",
            source="SDBS",
        )

    def get_spectrum_by_id(self, spectrum_id: str) -> Optional[ExperimentalSpectrum]:
        """Retrieve spectrum by SDBS ID (stub)."""
        logger.warning("SDBS spectrum retrieval not implemented (no public API)")
        return None

    def get_citation(self) -> str:
        """Get citation for SDBS."""
        return (
            "Spectral Database for Organic Compounds (SDBS). "
            "National Institute of Advanced Industrial Science and Technology (AIST), Japan. "
            "https://sdbs.db.aist.go.jp"
        )
