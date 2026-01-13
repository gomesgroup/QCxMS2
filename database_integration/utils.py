"""
Utility functions for mass spectrum processing and comparison.
"""

import math
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np

from .data_classes import ExperimentalSpectrum, Peak, SpectrumComparisonResult
from .enums import SimilarityMetric


def compare_spectra(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    metric: SimilarityMetric = SimilarityMetric.COSINE,
    mz_tolerance: float = 0.5,
    query_id: str = "query",
    reference_id: str = "reference",
) -> SpectrumComparisonResult:
    """
    Compare two mass spectra and calculate similarity.

    Args:
        query_peaks: Peaks from query spectrum
        reference_peaks: Peaks from reference spectrum
        metric: Similarity metric to use
        mz_tolerance: m/z tolerance for peak matching (Da)
        query_id: ID for query spectrum
        reference_id: ID for reference spectrum

    Returns:
        SpectrumComparisonResult with similarity score and matched peaks
    """
    # Match peaks within tolerance
    matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)

    # Calculate similarity based on metric
    if metric == SimilarityMetric.COSINE:
        score = cosine_similarity(query_peaks, reference_peaks, matched_pairs)
    elif metric == SimilarityMetric.DOT_PRODUCT:
        score = dot_product_similarity(query_peaks, reference_peaks, matched_pairs)
    elif metric == SimilarityMetric.EUCLIDEAN:
        score = euclidean_similarity(query_peaks, reference_peaks, matched_pairs)
    else:
        raise ValueError(f"Unsupported similarity metric: {metric}")

    return SpectrumComparisonResult(
        similarity_score=score,
        metric=metric,
        matched_peaks=len(matched_pairs),
        total_peaks_query=len(query_peaks),
        total_peaks_reference=len(reference_peaks),
        query_spectrum_id=query_id,
        reference_spectrum_id=reference_id,
        matched_peak_pairs=matched_pairs,
        mz_tolerance=mz_tolerance,
    )


def match_peaks(
    query_peaks: List[Peak], reference_peaks: List[Peak], tolerance: float = 0.5
) -> List[Tuple[Peak, Peak]]:
    """
    Match peaks between two spectra within a given m/z tolerance.

    Args:
        query_peaks: Peaks from query spectrum
        reference_peaks: Peaks from reference spectrum
        tolerance: m/z tolerance for matching (Da)

    Returns:
        List of (query_peak, reference_peak) tuples for matched peaks
    """
    matched_pairs = []

    # Sort peaks by m/z for efficient matching
    sorted_query = sorted(query_peaks, key=lambda p: p.mz)
    sorted_ref = sorted(reference_peaks, key=lambda p: p.mz)

    ref_idx = 0
    for q_peak in sorted_query:
        # Find matching reference peaks within tolerance
        best_match = None
        best_mz_diff = tolerance

        # Search forward from current position
        while ref_idx < len(sorted_ref):
            r_peak = sorted_ref[ref_idx]
            mz_diff = abs(q_peak.mz - r_peak.mz)

            if mz_diff <= tolerance and mz_diff < best_mz_diff:
                best_match = r_peak
                best_mz_diff = mz_diff

            # If reference peak is beyond tolerance, stop searching
            if r_peak.mz > q_peak.mz + tolerance:
                break

            ref_idx += 1

        if best_match is not None:
            matched_pairs.append((q_peak, best_match))
            # Reset index for next query peak
            ref_idx = max(0, ref_idx - 5)  # Backtrack a bit for overlapping matches

    return matched_pairs


def cosine_similarity(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    matched_pairs: List[Tuple[Peak, Peak]],
) -> float:
    """
    Calculate cosine similarity between two spectra.

    Cosine similarity = sum(q_i * r_i) / (sqrt(sum(q_i^2)) * sqrt(sum(r_i^2)))

    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        matched_pairs: Matched peak pairs

    Returns:
        Cosine similarity score (0-1)
    """
    if not matched_pairs:
        return 0.0

    # Calculate dot product of matched peaks
    dot_product = sum(q.intensity * r.intensity for q, r in matched_pairs)

    # Calculate magnitudes
    query_magnitude = math.sqrt(sum(p.intensity**2 for p in query_peaks))
    ref_magnitude = math.sqrt(sum(p.intensity**2 for p in reference_peaks))

    if query_magnitude == 0 or ref_magnitude == 0:
        return 0.0

    return dot_product / (query_magnitude * ref_magnitude)


def dot_product_similarity(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    matched_pairs: List[Tuple[Peak, Peak]],
) -> float:
    """
    Calculate simple dot product similarity (unnormalized).

    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        matched_pairs: Matched peak pairs

    Returns:
        Dot product score (normalized to 0-1 based on max possible)
    """
    if not matched_pairs:
        return 0.0

    dot_product = sum(q.intensity * r.intensity for q, r in matched_pairs)

    # Normalize by maximum possible dot product
    max_dot = min(
        sum(p.intensity**2 for p in query_peaks),
        sum(p.intensity**2 for p in reference_peaks),
    )

    if max_dot == 0:
        return 0.0

    return min(1.0, dot_product / max_dot)


def euclidean_similarity(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    matched_pairs: List[Tuple[Peak, Peak]],
) -> float:
    """
    Calculate Euclidean distance-based similarity (converted to 0-1 scale).

    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        matched_pairs: Matched peak pairs

    Returns:
        Similarity score (0-1, where 1 is identical)
    """
    if not matched_pairs:
        return 0.0

    # Calculate sum of squared differences for matched peaks
    squared_diff = sum((q.intensity - r.intensity) ** 2 for q, r in matched_pairs)

    # Add penalty for unmatched peaks
    matched_query_mzs = {q.mz for q, _ in matched_pairs}
    matched_ref_mzs = {r.mz for _, r in matched_pairs}

    unmatched_query = [p for p in query_peaks if p.mz not in matched_query_mzs]
    unmatched_ref = [p for p in reference_peaks if p.mz not in matched_ref_mzs]

    squared_diff += sum(p.intensity**2 for p in unmatched_query)
    squared_diff += sum(p.intensity**2 for p in unmatched_ref)

    distance = math.sqrt(squared_diff)

    # Convert distance to similarity (0-1 scale)
    # Max distance is sqrt(sum of all intensities squared)
    max_intensity = max(
        math.sqrt(sum(p.intensity**2 for p in query_peaks)),
        math.sqrt(sum(p.intensity**2 for p in reference_peaks)),
    )

    if max_intensity == 0:
        return 0.0

    return max(0.0, 1.0 - (distance / (max_intensity * math.sqrt(2))))


def parse_spectrum_file(
    file_path: Union[str, Path], format_type: str = "auto"
) -> ExperimentalSpectrum:
    """
    Parse a mass spectrum file into ExperimentalSpectrum.

    Supported formats:
    - MSP (NIST format)
    - MGF (Mascot Generic Format)
    - TXT (simple two-column format: mz intensity)

    Args:
        file_path: Path to spectrum file
        format_type: Format type ("msp", "mgf", "txt", or "auto" for auto-detection)

    Returns:
        ExperimentalSpectrum object

    Raises:
        ValueError: If format is not supported or file cannot be parsed
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise ValueError(f"File not found: {file_path}")

    # Auto-detect format from extension
    if format_type == "auto":
        ext = file_path.suffix.lower()
        if ext == ".msp":
            format_type = "msp"
        elif ext == ".mgf":
            format_type = "mgf"
        elif ext in [".txt", ".dat"]:
            format_type = "txt"
        else:
            raise ValueError(f"Unknown file format: {ext}. Specify format_type explicitly.")

    # Parse based on format
    if format_type == "msp":
        return _parse_msp_file(file_path)
    elif format_type == "mgf":
        return _parse_mgf_file(file_path)
    elif format_type == "txt":
        return _parse_txt_file(file_path)
    else:
        raise ValueError(f"Unsupported format: {format_type}")


def _parse_txt_file(file_path: Path) -> ExperimentalSpectrum:
    """Parse simple two-column text file (mz intensity)."""
    peaks = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    peaks.append(Peak(mz=mz, intensity=intensity))
                except ValueError:
                    continue

    return ExperimentalSpectrum(
        spectrum_id=file_path.stem,
        name=file_path.stem,
        peaks=peaks,
        source="local_file",
    )


def _parse_msp_file(file_path: Path) -> ExperimentalSpectrum:
    """Parse MSP format file (NIST format) - basic implementation."""
    # Simplified parser - full MSP parsing is more complex
    with open(file_path, "r") as f:
        lines = f.readlines()

    metadata = {}
    peaks = []
    in_peaks = False

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if ":" in line and not in_peaks:
            key, value = line.split(":", 1)
            metadata[key.strip()] = value.strip()
        elif line.startswith("Num peaks"):
            in_peaks = True
        elif in_peaks:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    peaks.append(Peak(mz=mz, intensity=intensity))
                except ValueError:
                    continue

    return ExperimentalSpectrum(
        spectrum_id=file_path.stem,
        name=metadata.get("Name", file_path.stem),
        formula=metadata.get("Formula"),
        peaks=peaks,
        source="msp_file",
    )


def _parse_mgf_file(file_path: Path) -> ExperimentalSpectrum:
    """Parse MGF format file - basic implementation."""
    # Simplified MGF parser
    with open(file_path, "r") as f:
        lines = f.readlines()

    metadata = {}
    peaks = []
    in_spectrum = False

    for line in lines:
        line = line.strip()

        if line == "BEGIN IONS":
            in_spectrum = True
            continue
        elif line == "END IONS":
            break
        elif not in_spectrum:
            continue

        if "=" in line:
            key, value = line.split("=", 1)
            metadata[key.strip()] = value.strip()
        else:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    peaks.append(Peak(mz=mz, intensity=intensity))
                except ValueError:
                    continue

    return ExperimentalSpectrum(
        spectrum_id=metadata.get("TITLE", file_path.stem),
        name=metadata.get("TITLE", file_path.stem),
        precursor_mz=float(metadata.get("PEPMASS", 0)),
        peaks=peaks,
        source="mgf_file",
    )
