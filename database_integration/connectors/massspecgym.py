"""
MassSpecGym connector for accessing the NeurIPS 2024 Spotlight MS/MS benchmark dataset.

MassSpecGym is a curated dataset of 231,104 high-quality MS/MS spectra from 29K molecules,
combining data from MassBank, MoNA, and GNPS with pre-defined train/val/test splits.

Paper: https://arxiv.org/abs/2410.23326
GitHub: https://github.com/pluskal-lab/MassSpecGym
Hugging Face: https://huggingface.co/datasets/roman-bushuiev/MassSpecGym
"""

import hashlib
import json
import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Union

from ..base import MSDatabaseConnector
from ..config import DatabaseConfig
from ..data_classes import ExperimentalSpectrum, MSSearchResult, Peak
from ..enums import DatabaseType, IonizationMode, SpectrumType
from ..search_criteria import SpectrumSearchCriteria

logger = logging.getLogger(__name__)

# Dataset configuration
DATASET_NAME = "roman-bushuiev/MassSpecGym"
DATASET_SPLITS = ["train", "val", "test"]


class MassSpecGymConnector(MSDatabaseConnector):
    """
    Connector for MassSpecGym benchmark dataset (NeurIPS 2024 Spotlight).

    Features:
    - 231,104 curated MS/MS spectra from 29K molecules
    - Pre-defined train/val/test splits for ML benchmarking
    - Combines MassBank, MoNA, GNPS into one dataset
    - Local caching for fast repeated access
    - Search by formula, InChI key, SMILES, or spectrum similarity

    Note: This connector uses the Hugging Face datasets library which downloads
    the dataset on first access (~2GB). Subsequent accesses use the local cache.
    """

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize MassSpecGym connector.

        Args:
            config: Optional custom configuration
        """
        if config is None:
            config = DatabaseConfig.for_massspecgym()
        super().__init__(config)

        self._dataset = None
        self._index_by_formula: Dict[str, List[int]] = {}
        self._index_by_inchikey: Dict[str, List[int]] = {}
        self._index_by_smiles: Dict[str, List[int]] = {}
        self._index_built = False
        self._splits_to_use = config.additional_params.get("splits", DATASET_SPLITS)

    def _ensure_dataset_loaded(self) -> bool:
        """
        Ensure the dataset is loaded. Downloads on first access.

        Returns:
            True if dataset is available, False otherwise
        """
        if self._dataset is not None:
            return True

        try:
            from datasets import load_dataset, concatenate_datasets

            logger.info(f"Loading MassSpecGym dataset from Hugging Face...")
            logger.info(f"This may take a few minutes on first download (~2GB)")

            # Load specified splits
            datasets_list = []
            for split in self._splits_to_use:
                try:
                    ds = load_dataset(DATASET_NAME, split=split)
                    datasets_list.append(ds)
                    logger.info(f"  Loaded {split} split: {len(ds)} spectra")
                except Exception as e:
                    logger.warning(f"  Failed to load {split} split: {e}")

            if not datasets_list:
                logger.error("No splits could be loaded")
                return False

            # Concatenate all splits
            if len(datasets_list) == 1:
                self._dataset = datasets_list[0]
            else:
                self._dataset = concatenate_datasets(datasets_list)

            logger.info(f"Total spectra loaded: {len(self._dataset)}")
            return True

        except ImportError:
            logger.error(
                "The 'datasets' package is required for MassSpecGym. "
                "Install it with: pip install datasets"
            )
            return False
        except Exception as e:
            logger.error(f"Failed to load MassSpecGym dataset: {e}")
            return False

    def _build_indexes(self) -> None:
        """Build indexes for fast searching by formula, InChI key, and SMILES."""
        if self._index_built or self._dataset is None:
            return

        logger.info("Building search indexes for MassSpecGym...")
        start_time = time.time()

        for idx in range(len(self._dataset)):
            try:
                record = self._dataset[idx]

                # Index by formula
                formula = record.get("formula", "")
                if formula:
                    if formula not in self._index_by_formula:
                        self._index_by_formula[formula] = []
                    self._index_by_formula[formula].append(idx)

                # Index by InChI key (2D key)
                inchikey = record.get("inchikey", "")
                if inchikey:
                    # Use first 14 chars (connectivity layer) for broader matching
                    inchikey_2d = inchikey[:14] if len(inchikey) >= 14 else inchikey
                    if inchikey_2d not in self._index_by_inchikey:
                        self._index_by_inchikey[inchikey_2d] = []
                    self._index_by_inchikey[inchikey_2d].append(idx)
                    # Also index full key
                    if inchikey not in self._index_by_inchikey:
                        self._index_by_inchikey[inchikey] = []
                    self._index_by_inchikey[inchikey].append(idx)

                # Index by SMILES (canonicalized if possible)
                smiles = record.get("smiles", "")
                if smiles:
                    if smiles not in self._index_by_smiles:
                        self._index_by_smiles[smiles] = []
                    self._index_by_smiles[smiles].append(idx)

            except Exception as e:
                logger.debug(f"Error indexing record {idx}: {e}")
                continue

        self._index_built = True
        elapsed = time.time() - start_time
        logger.info(
            f"Indexes built in {elapsed:.1f}s: "
            f"{len(self._index_by_formula)} formulas, "
            f"{len(self._index_by_inchikey)} InChI keys, "
            f"{len(self._index_by_smiles)} SMILES"
        )

    def test_connection(self) -> bool:
        """Test connection to MassSpecGym (verify dataset can be loaded)."""
        try:
            return self._ensure_dataset_loaded()
        except Exception as e:
            logger.error(f"MassSpecGym connection test failed: {e}")
            return False

    def search_spectrum(self, criteria: SpectrumSearchCriteria) -> MSSearchResult:
        """
        Search MassSpecGym for spectra matching criteria.

        Args:
            criteria: Search criteria

        Returns:
            MSSearchResult with matching spectra
        """
        criteria.validate()

        if not self._ensure_dataset_loaded():
            return MSSearchResult(
                success=False,
                error_message="Failed to load MassSpecGym dataset. Install 'datasets' package.",
                source="MassSpecGym",
            )

        # Build indexes if not already done
        self._build_indexes()

        start_time = time.time()

        # Determine search method based on criteria
        if criteria.inchi_key:
            indices = self._search_by_inchikey(criteria.inchi_key)
        elif criteria.smiles:
            indices = self._search_by_smiles(criteria.smiles)
        elif criteria.formula:
            indices = self._search_by_formula(criteria.formula)
        elif criteria.peaks:
            indices = self._search_by_peaks(criteria)
        else:
            return MSSearchResult(
                success=False,
                error_message="MassSpecGym requires formula, InChI key, SMILES, or peaks for search",
                source="MassSpecGym",
            )

        # Apply ionization mode filter if specified
        if criteria.ionization_mode:
            indices = self._filter_by_ionization(indices, criteria.ionization_mode)

        # Apply max results limit
        if criteria.max_results and len(indices) > criteria.max_results:
            indices = indices[: criteria.max_results]

        # Convert to ExperimentalSpectrum objects
        spectra = []
        for idx in indices:
            try:
                spectrum = self._parse_record(self._dataset[idx], idx)
                if spectrum:
                    spectra.append(spectrum)
            except Exception as e:
                logger.debug(f"Error parsing record {idx}: {e}")
                continue

        request_time = time.time() - start_time

        return MSSearchResult(
            success=True,
            spectra=spectra,
            total_results=len(spectra),
            source="MassSpecGym",
            citation=self.get_citation(),
            request_time=request_time,
            cached=False,  # Dataset is cached locally by HF
        )

    def _search_by_formula(self, formula: str) -> List[int]:
        """Search by molecular formula."""
        return self._index_by_formula.get(formula, [])

    def _search_by_inchikey(self, inchikey: str) -> List[int]:
        """Search by InChI key (supports full key or 2D connectivity layer)."""
        # Try full key first
        if inchikey in self._index_by_inchikey:
            return self._index_by_inchikey[inchikey]

        # Try 2D key (first 14 chars)
        inchikey_2d = inchikey[:14] if len(inchikey) >= 14 else inchikey
        return self._index_by_inchikey.get(inchikey_2d, [])

    def _search_by_smiles(self, smiles: str) -> List[int]:
        """Search by SMILES string."""
        return self._index_by_smiles.get(smiles, [])

    def _search_by_peaks(self, criteria: SpectrumSearchCriteria) -> List[int]:
        """
        Search by peak list (similarity search).

        This is a computationally expensive operation for large datasets.
        Consider using pre-built similarity indexes for production use.
        """
        if not criteria.peaks:
            return []

        query_peaks = [Peak(mz=mz, intensity=intensity) for mz, intensity in criteria.peaks]
        mz_tolerance = criteria.mz_tolerance
        threshold = criteria.similarity_threshold

        # Limit search to avoid excessive computation
        max_scan = min(10000, len(self._dataset))
        logger.info(f"Performing peak similarity search on first {max_scan} spectra...")

        matches = []
        for idx in range(max_scan):
            try:
                record = self._dataset[idx]
                ref_mzs = record.get("mzs", [])
                ref_intensities = record.get("intensities", [])

                if not ref_mzs or not ref_intensities:
                    continue

                ref_peaks = [Peak(mz=mz, intensity=intensity) for mz, intensity in zip(ref_mzs, ref_intensities)]

                # Calculate cosine similarity
                similarity = self._calculate_similarity(query_peaks, ref_peaks, mz_tolerance)

                if similarity >= threshold:
                    matches.append((idx, similarity))

            except Exception as e:
                logger.debug(f"Error comparing spectrum {idx}: {e}")
                continue

        # Sort by similarity
        matches.sort(key=lambda x: x[1], reverse=True)
        return [idx for idx, _ in matches]

    def _calculate_similarity(
        self, query_peaks: List[Peak], ref_peaks: List[Peak], mz_tolerance: float
    ) -> float:
        """Calculate cosine similarity between two peak lists."""
        if not query_peaks or not ref_peaks:
            return 0.0

        # Match peaks within tolerance
        matched_query = []
        matched_ref = []

        for qp in query_peaks:
            best_match = None
            best_diff = mz_tolerance + 1

            for rp in ref_peaks:
                diff = abs(qp.mz - rp.mz)
                if diff <= mz_tolerance and diff < best_diff:
                    best_match = rp
                    best_diff = diff

            if best_match:
                matched_query.append(qp.intensity)
                matched_ref.append(best_match.intensity)

        if not matched_query:
            return 0.0

        # Cosine similarity
        import math

        dot_product = sum(q * r for q, r in zip(matched_query, matched_ref))
        norm_q = math.sqrt(sum(q**2 for q in matched_query))
        norm_r = math.sqrt(sum(r**2 for r in matched_ref))

        if norm_q == 0 or norm_r == 0:
            return 0.0

        return dot_product / (norm_q * norm_r)

    def _filter_by_ionization(self, indices: List[int], mode: IonizationMode) -> List[int]:
        """Filter indices by ionization mode."""
        filtered = []
        for idx in indices:
            try:
                record = self._dataset[idx]
                adduct = record.get("adduct", "")

                # Map adduct to ionization mode
                record_mode = self._parse_adduct_to_mode(adduct)
                if record_mode == mode or record_mode == IonizationMode.OTHER:
                    filtered.append(idx)
            except Exception:
                continue
        return filtered

    @staticmethod
    def _parse_adduct_to_mode(adduct: str) -> IonizationMode:
        """Parse adduct string to ionization mode."""
        if not adduct:
            return IonizationMode.OTHER

        adduct_upper = adduct.upper()
        if "+H]+" in adduct_upper or "+NA]+" in adduct_upper or "+K]+" in adduct_upper:
            return IonizationMode.ESI_POSITIVE
        elif "-H]-" in adduct_upper or "+CL]-" in adduct_upper:
            return IonizationMode.ESI_NEGATIVE
        elif "[M]+" in adduct_upper or "[M]" in adduct_upper:
            return IonizationMode.EI  # Often EI for radical cation

        return IonizationMode.OTHER

    def get_spectrum_by_id(self, spectrum_id: str) -> Optional[ExperimentalSpectrum]:
        """
        Retrieve spectrum by index or identifier.

        Args:
            spectrum_id: Index in dataset (as string) or InChI key

        Returns:
            ExperimentalSpectrum or None if not found
        """
        if not self._ensure_dataset_loaded():
            return None

        # Try as integer index
        try:
            idx = int(spectrum_id)
            if 0 <= idx < len(self._dataset):
                return self._parse_record(self._dataset[idx], idx)
        except ValueError:
            pass

        # Try as InChI key
        self._build_indexes()
        indices = self._search_by_inchikey(spectrum_id)
        if indices:
            return self._parse_record(self._dataset[indices[0]], indices[0])

        return None

    def _parse_record(self, record: Dict[str, Any], idx: int) -> Optional[ExperimentalSpectrum]:
        """Parse a MassSpecGym record into ExperimentalSpectrum."""
        try:
            # Extract peaks
            mzs = record.get("mzs", [])
            intensities = record.get("intensities", [])

            if not mzs or not intensities:
                return None

            peaks = [Peak(mz=float(mz), intensity=float(intensity)) for mz, intensity in zip(mzs, intensities)]

            # Extract metadata
            smiles = record.get("smiles", "")
            inchikey = record.get("inchikey", "")
            formula = record.get("formula", "")
            precursor_mz = record.get("precursor_mz", None)
            adduct = record.get("adduct", "")
            instrument_type = record.get("instrument_type", "")
            collision_energy = record.get("collision_energy", "")

            # Generate a stable spectrum ID from InChI key and index
            spectrum_id = f"MSG-{inchikey[:14] if inchikey else 'UNK'}-{idx}"

            # Parse ionization mode from adduct
            ionization_mode = self._parse_adduct_to_mode(adduct)

            spectrum = ExperimentalSpectrum(
                spectrum_id=spectrum_id,
                name=smiles[:50] if smiles else f"MassSpecGym-{idx}",
                formula=formula,
                smiles=smiles,
                inchi_key=inchikey,
                peaks=peaks,
                precursor_mz=float(precursor_mz) if precursor_mz else None,
                precursor_type=adduct,
                ionization_mode=ionization_mode,
                spectrum_type=SpectrumType.MS2,
                collision_energy=str(collision_energy) if collision_energy else None,
                instrument_type=instrument_type,
                source="MassSpecGym",
                database_id=str(idx),
                url=f"https://huggingface.co/datasets/roman-bushuiev/MassSpecGym",
                metadata={
                    "dataset_index": idx,
                    "adduct": adduct,
                    "original_source": record.get("source", ""),
                },
            )

            return spectrum

        except Exception as e:
            logger.error(f"Error parsing MassSpecGym record {idx}: {e}")
            return None

    def get_dataset_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about the loaded dataset.

        Returns:
            Dictionary with dataset statistics
        """
        if not self._ensure_dataset_loaded():
            return {"error": "Dataset not loaded"}

        self._build_indexes()

        return {
            "total_spectra": len(self._dataset),
            "unique_formulas": len(self._index_by_formula),
            "unique_inchikeys": len(self._index_by_inchikey),
            "unique_smiles": len(self._index_by_smiles),
            "splits_loaded": self._splits_to_use,
            "source": "MassSpecGym (NeurIPS 2024 Spotlight)",
        }

    def iterate_spectra(
        self, max_spectra: Optional[int] = None, start_idx: int = 0
    ) -> Iterator[ExperimentalSpectrum]:
        """
        Iterate over spectra in the dataset.

        Args:
            max_spectra: Maximum number of spectra to yield (None for all)
            start_idx: Starting index

        Yields:
            ExperimentalSpectrum objects
        """
        if not self._ensure_dataset_loaded():
            return

        end_idx = len(self._dataset)
        if max_spectra:
            end_idx = min(start_idx + max_spectra, len(self._dataset))

        for idx in range(start_idx, end_idx):
            try:
                spectrum = self._parse_record(self._dataset[idx], idx)
                if spectrum:
                    yield spectrum
            except Exception as e:
                logger.debug(f"Error iterating spectrum {idx}: {e}")
                continue

    def get_citation(self) -> str:
        """Get citation for MassSpecGym."""
        return (
            "Bushuiev et al. MassSpecGym: A benchmark for the discovery and identification "
            "of molecules. NeurIPS 2024 (Spotlight). "
            "https://github.com/pluskal-lab/MassSpecGym"
        )

    # MSDatabaseConnector does not require _make_request for local dataset
    def _make_request(self, url: str, **kwargs) -> Dict[str, Any]:
        """Not used - MassSpecGym uses local Hugging Face dataset."""
        raise NotImplementedError("MassSpecGym uses local dataset, not HTTP API")
