"""
Base connector class for mass spectrometry databases.

Simplified version inspired by Explorer's database integration architecture.
"""

import json
import logging
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import urlencode, urlparse

try:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
except ImportError:
    raise ImportError(
        "requests library is required. Install with: pip install requests"
    )

from .config import DatabaseConfig
from .data_classes import ExperimentalSpectrum, MSSearchResult
from .search_criteria import SpectrumSearchCriteria

logger = logging.getLogger(__name__)


class MSDatabaseConnector(ABC):
    """
    Abstract base class for mass spectrometry database connectors.

    Provides common functionality for HTTP requests, caching, rate limiting,
    and error handling. Database-specific connectors inherit from this class.
    """

    def __init__(self, config: DatabaseConfig):
        """
        Initialize connector with configuration.

        Args:
            config: Database configuration
        """
        self.config = config
        self.session = self._create_session()
        self._last_request_time = 0.0
        self._request_count = 0

    def _create_session(self) -> requests.Session:
        """
        Create and configure HTTP session with retries.

        Returns:
            Configured requests Session
        """
        session = requests.Session()

        # Configure retries
        retry_strategy = Retry(
            total=self.config.max_retries,
            backoff_factor=self.config.retry_delay,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["HEAD", "GET", "OPTIONS", "POST"],
        )

        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)

        # Set user agent
        session.headers.update({"User-Agent": self.config.user_agent})

        # Add API key if provided
        if self.config.api_key:
            session.headers.update({"X-API-Key": self.config.api_key})

        # Add basic auth if provided
        if self.config.username and self.config.password:
            session.auth = (self.config.username, self.config.password)

        return session

    def _rate_limit(self):
        """Apply rate limiting between requests."""
        if self.config.rate_limit_per_second <= 0:
            return

        min_interval = 1.0 / self.config.rate_limit_per_second
        elapsed = time.time() - self._last_request_time

        if elapsed < min_interval:
            sleep_time = min_interval - elapsed
            logger.debug(f"Rate limiting: sleeping for {sleep_time:.2f}s")
            time.sleep(sleep_time)

    def _make_request(
        self,
        url: str,
        method: str = "GET",
        params: Optional[Dict[str, Any]] = None,
        data: Optional[Dict[str, Any]] = None,
        json_data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """
        Make HTTP request with rate limiting, caching, and error handling.

        Args:
            url: URL to request
            method: HTTP method (GET, POST, etc.)
            params: URL query parameters
            data: Form data for POST
            json_data: JSON data for POST
            headers: Additional headers
            use_cache: Whether to use cached results

        Returns:
            Dictionary with success, data, error_message, etc.
        """
        # Check cache first
        if use_cache and self.config.cache_enabled and method == "GET":
            cache_key = self._get_cache_key(url, params)
            cached_result = self._load_from_cache(cache_key)
            if cached_result is not None:
                logger.debug(f"Cache hit for {url}")
                return {
                    "success": True,
                    "data": cached_result,
                    "cached": True,
                    "status_code": 200,
                    "request_time": 0.0,
                }

        # Apply rate limiting
        self._rate_limit()

        # Make request
        start_time = time.time()
        try:
            response = self.session.request(
                method=method,
                url=url,
                params=params,
                data=data,
                json=json_data,
                headers=headers,
                timeout=self.config.timeout,
            )

            self._last_request_time = time.time()
            self._request_count += 1
            request_time = time.time() - start_time

            # Handle response
            response.raise_for_status()

            # Try to parse as JSON, fallback to text
            try:
                data_result = response.json()
            except json.JSONDecodeError:
                data_result = response.text

            # Cache successful GET requests
            if use_cache and self.config.cache_enabled and method == "GET":
                cache_key = self._get_cache_key(url, params)
                self._save_to_cache(cache_key, data_result)

            return {
                "success": True,
                "data": data_result,
                "cached": False,
                "status_code": response.status_code,
                "request_time": request_time,
                "headers": dict(response.headers),
            }

        except requests.exceptions.HTTPError as e:
            logger.error(f"HTTP error for {url}: {e}")
            return {
                "success": False,
                "error_message": f"HTTP {e.response.status_code}: {e.response.reason}",
                "status_code": e.response.status_code,
                "request_time": time.time() - start_time,
            }

        except requests.exceptions.Timeout:
            logger.error(f"Timeout for {url}")
            return {
                "success": False,
                "error_message": f"Request timeout after {self.config.timeout}s",
                "status_code": 408,
                "request_time": time.time() - start_time,
            }

        except requests.exceptions.RequestException as e:
            logger.error(f"Request error for {url}: {e}")
            return {
                "success": False,
                "error_message": str(e),
                "status_code": 0,
                "request_time": time.time() - start_time,
            }

    def _get_cache_key(self, url: str, params: Optional[Dict[str, Any]]) -> str:
        """Generate cache key from URL and parameters."""
        if params:
            url_with_params = f"{url}?{urlencode(sorted(params.items()))}"
        else:
            url_with_params = url

        # Use hash of URL as filename
        import hashlib

        return hashlib.md5(url_with_params.encode()).hexdigest()

    def _load_from_cache(self, cache_key: str) -> Optional[Any]:
        """Load data from cache if not expired."""
        if not self.config.cache_enabled:
            return None

        cache_file = Path(self.config.cache_dir) / f"{cache_key}.json"

        if not cache_file.exists():
            return None

        # Check if cache is expired
        age = time.time() - cache_file.stat().st_mtime
        if age > self.config.cache_ttl:
            logger.debug(f"Cache expired for {cache_key}")
            cache_file.unlink()
            return None

        try:
            with open(cache_file, "r") as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Failed to load cache {cache_key}: {e}")
            return None

    def _save_to_cache(self, cache_key: str, data: Any):
        """Save data to cache."""
        if not self.config.cache_enabled:
            return

        cache_file = Path(self.config.cache_dir) / f"{cache_key}.json"

        try:
            with open(cache_file, "w") as f:
                json.dump(data, f)
        except (TypeError, IOError) as e:
            logger.warning(f"Failed to save cache {cache_key}: {e}")

    @abstractmethod
    def test_connection(self) -> bool:
        """
        Test connection to the database.

        Returns:
            True if connection successful, False otherwise
        """
        pass

    @abstractmethod
    def search_spectrum(
        self, criteria: SpectrumSearchCriteria
    ) -> MSSearchResult:
        """
        Search for mass spectra matching criteria.

        Args:
            criteria: Search criteria

        Returns:
            MSSearchResult with matching spectra
        """
        pass

    @abstractmethod
    def get_spectrum_by_id(self, spectrum_id: str) -> Optional[ExperimentalSpectrum]:
        """
        Retrieve a specific spectrum by database ID.

        Args:
            spectrum_id: Database-specific spectrum ID

        Returns:
            ExperimentalSpectrum or None if not found
        """
        pass

    def get_citation(self) -> str:
        """
        Get citation information for the database.

        Returns:
            Citation string
        """
        return f"Data from {self.config.name}"

    def close(self):
        """Close the HTTP session."""
        if self.session:
            self.session.close()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
