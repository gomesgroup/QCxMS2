"""
Configuration for mass spectrometry database connectors.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional

from .enums import DatabaseType

# Default configuration values
DEFAULT_TIMEOUT = 30  # seconds
DEFAULT_MAX_RETRIES = 3
DEFAULT_RETRY_DELAY = 1.0  # seconds
DEFAULT_CACHE_DIR = os.path.expanduser("~/.cache/qcxms2/ms_databases")
DEFAULT_USER_AGENT = "QCxMS2-Database-Integration/1.0"


@dataclass
class DatabaseConfig:
    """
    Configuration for a mass spectrometry database connector.

    Attributes:
        database_type: Type of database (MassBank, NIST, SDBS, etc.)
        name: Human-readable name for the database
        base_url: Base URL for API/web access
        api_key: Optional API key for authenticated access
        username: Optional username for authentication
        password: Optional password for authentication
        timeout: Request timeout in seconds
        max_retries: Maximum number of retry attempts
        retry_delay: Delay between retries in seconds
        cache_dir: Directory for caching results
        cache_enabled: Whether to enable caching
        cache_ttl: Cache time-to-live in seconds
        user_agent: User agent string for HTTP requests
        rate_limit_per_second: Maximum requests per second
        additional_params: Additional database-specific parameters
    """

    database_type: DatabaseType
    name: str
    base_url: Optional[str] = None
    api_key: Optional[str] = None
    username: Optional[str] = None
    password: Optional[str] = None
    timeout: int = DEFAULT_TIMEOUT
    max_retries: int = DEFAULT_MAX_RETRIES
    retry_delay: float = DEFAULT_RETRY_DELAY
    cache_dir: str = DEFAULT_CACHE_DIR
    cache_enabled: bool = True
    cache_ttl: int = 86400  # 24 hours
    user_agent: str = DEFAULT_USER_AGENT
    rate_limit_per_second: float = 2.0  # Conservative default
    additional_params: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Initialize cache directory."""
        if self.cache_enabled:
            cache_path = Path(self.cache_dir) / self.database_type.value / self.name
            cache_path.mkdir(parents=True, exist_ok=True)
            self.cache_dir = str(cache_path)

    @classmethod
    def for_massbank(cls, **kwargs) -> "DatabaseConfig":
        """Create configuration for MassBank."""
        defaults = {
            "database_type": DatabaseType.MASSBANK,
            "name": "massbank_eu",
            "base_url": "https://massbank.eu/MassBank/api",
            "rate_limit_per_second": 2.0,  # Be conservative
        }
        defaults.update(kwargs)
        return cls(**defaults)

    @classmethod
    def for_mona(cls, **kwargs) -> "DatabaseConfig":
        """Create configuration for MassBank of North America (MoNA)."""
        defaults = {
            "database_type": DatabaseType.MONA,
            "name": "mona",
            "base_url": "https://mona.fiehnlab.ucdavis.edu/rest",
            "rate_limit_per_second": 2.0,
        }
        defaults.update(kwargs)
        return cls(**defaults)

    @classmethod
    def for_nist_ms(cls, **kwargs) -> "DatabaseConfig":
        """Create configuration for NIST Chemistry WebBook."""
        defaults = {
            "database_type": DatabaseType.NIST_MS,
            "name": "nist_webbook",
            "base_url": "https://webbook.nist.gov/cgi",
            "rate_limit_per_second": 1.0,  # NIST might be more restrictive
        }
        defaults.update(kwargs)
        return cls(**defaults)

    @classmethod
    def for_sdbs(cls, **kwargs) -> "DatabaseConfig":
        """Create configuration for SDBS (Spectral Database for Organic Compounds)."""
        defaults = {
            "database_type": DatabaseType.SDBS,
            "name": "sdbs_aist",
            "base_url": "https://sdbs.db.aist.go.jp/sdbs/cgi-bin",
            "rate_limit_per_second": 0.5,  # SDBS has 50/day limit, be very conservative
            "additional_params": {"daily_limit": 50},
        }
        defaults.update(kwargs)
        return cls(**defaults)

    @classmethod
    def for_gnps(cls, **kwargs) -> "DatabaseConfig":
        """Create configuration for GNPS."""
        defaults = {
            "database_type": DatabaseType.GNPS,
            "name": "gnps",
            "base_url": "https://gnps.ucsd.edu",
            "rate_limit_per_second": 2.0,
        }
        defaults.update(kwargs)
        return cls(**defaults)

    @classmethod
    def for_massspecgym(cls, **kwargs) -> "DatabaseConfig":
        """
        Create configuration for MassSpecGym (NeurIPS 2024 benchmark).

        MassSpecGym is a local dataset loaded via Hugging Face datasets library,
        so it doesn't have a base_url or rate limiting.
        """
        defaults = {
            "database_type": DatabaseType.MASSSPECGYM,
            "name": "massspecgym",
            "base_url": None,  # Local dataset, no API
            "rate_limit_per_second": 0,  # No rate limiting needed
            "cache_enabled": True,  # HF datasets handles caching
            "additional_params": {
                "splits": ["train", "val", "test"],  # Which splits to load
                "hf_dataset": "roman-bushuiev/MassSpecGym",
            },
        }
        defaults.update(kwargs)
        return cls(**defaults)
