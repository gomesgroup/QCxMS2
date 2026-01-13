"""
Similarity metrics and statistical analysis for QCxMS2 benchmark.

This module provides:
- Various spectrum similarity metrics (cosine, dot product, spectral entropy)
- Statistical analysis functions (mean, std, percentiles, confidence intervals)
- Matching quality assessment (peak coverage, intensity correlation)

References:
- Stein & Scott (1994) - Cosine similarity for MS matching
- Li et al. (2021) - Spectral entropy for MS/MS similarity
- MassSpecGym (NeurIPS 2024) - Benchmark metrics
"""

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np


@dataclass
class Peak:
    """A single peak in a mass spectrum."""
    mz: float
    intensity: float
    annotation: Optional[str] = None


@dataclass
class SimilarityMetrics:
    """Container for multiple similarity metrics."""
    
    cosine: float = 0.0
    dot_product: float = 0.0
    spectral_entropy: float = 0.0
    weighted_cosine: float = 0.0
    reverse_cosine: float = 0.0
    matched_peaks_ratio: float = 0.0
    intensity_correlation: float = 0.0
    
    # Peak matching details
    matched_peaks: int = 0
    total_query_peaks: int = 0
    total_reference_peaks: int = 0
    
    def to_dict(self) -> Dict[str, float]:
        """Convert to dictionary."""
        return {
            "cosine": self.cosine,
            "dot_product": self.dot_product,
            "spectral_entropy": self.spectral_entropy,
            "weighted_cosine": self.weighted_cosine,
            "reverse_cosine": self.reverse_cosine,
            "matched_peaks_ratio": self.matched_peaks_ratio,
            "intensity_correlation": self.intensity_correlation,
            "matched_peaks": self.matched_peaks,
            "total_query_peaks": self.total_query_peaks,
            "total_reference_peaks": self.total_reference_peaks,
        }
    
    @property
    def average(self) -> float:
        """Calculate average of main metrics."""
        metrics = [self.cosine, self.dot_product, self.spectral_entropy]
        return sum(metrics) / len(metrics) if metrics else 0.0


@dataclass
class StatisticalAnalysis:
    """Statistical analysis results for a set of similarity scores."""
    
    count: int = 0
    mean: float = 0.0
    std: float = 0.0
    median: float = 0.0
    min_value: float = 0.0
    max_value: float = 0.0
    percentile_25: float = 0.0
    percentile_75: float = 0.0
    percentile_90: float = 0.0
    percentile_95: float = 0.0
    
    # Quality metrics
    above_threshold_50: float = 0.0  # % above 0.5
    above_threshold_70: float = 0.0  # % above 0.7
    above_threshold_90: float = 0.0  # % above 0.9
    
    # Confidence interval (95%)
    ci_lower: float = 0.0
    ci_upper: float = 0.0
    
    def to_dict(self) -> Dict[str, float]:
        """Convert to dictionary."""
        return {
            "count": self.count,
            "mean": self.mean,
            "std": self.std,
            "median": self.median,
            "min": self.min_value,
            "max": self.max_value,
            "percentile_25": self.percentile_25,
            "percentile_75": self.percentile_75,
            "percentile_90": self.percentile_90,
            "percentile_95": self.percentile_95,
            "above_threshold_50_pct": self.above_threshold_50,
            "above_threshold_70_pct": self.above_threshold_70,
            "above_threshold_90_pct": self.above_threshold_90,
            "ci_95_lower": self.ci_lower,
            "ci_95_upper": self.ci_upper,
        }
    
    def __str__(self) -> str:
        """Human-readable summary."""
        return (
            f"n={self.count}, mean={self.mean:.4f}±{self.std:.4f}, "
            f"median={self.median:.4f}, range=[{self.min_value:.4f}, {self.max_value:.4f}]"
        )


def match_peaks(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
) -> List[Tuple[Peak, Peak, float]]:
    """
    Match peaks between query and reference spectra within m/z tolerance.
    
    Uses a greedy algorithm that matches each query peak to the closest
    reference peak within tolerance, prioritizing higher intensity matches.
    
    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        mz_tolerance: m/z tolerance for matching (Da)
    
    Returns:
        List of (query_peak, reference_peak, mz_difference) tuples
    """
    if not query_peaks or not reference_peaks:
        return []
    
    # Sort peaks by m/z for efficient matching
    sorted_query = sorted(query_peaks, key=lambda p: p.mz)
    sorted_ref = sorted(reference_peaks, key=lambda p: p.mz)
    
    matched_pairs = []
    used_ref_indices = set()
    
    # Match query peaks to closest reference peaks
    for q_peak in sorted_query:
        best_match_idx = None
        best_mz_diff = mz_tolerance + 1
        
        for i, r_peak in enumerate(sorted_ref):
            if i in used_ref_indices:
                continue
            
            mz_diff = abs(q_peak.mz - r_peak.mz)
            
            if mz_diff <= mz_tolerance and mz_diff < best_mz_diff:
                best_match_idx = i
                best_mz_diff = mz_diff
            
            # Early termination if reference peak is beyond tolerance
            if r_peak.mz > q_peak.mz + mz_tolerance:
                break
        
        if best_match_idx is not None:
            matched_pairs.append((q_peak, sorted_ref[best_match_idx], best_mz_diff))
            used_ref_indices.add(best_match_idx)
    
    return matched_pairs


def calculate_cosine_similarity(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
    matched_pairs: Optional[List[Tuple[Peak, Peak, float]]] = None,
) -> float:
    """
    Calculate cosine similarity between two spectra.
    
    Cosine similarity is the standard metric for mass spectrum matching,
    measuring the angle between two intensity vectors.
    
    Formula: cos(θ) = Σ(qi × ri) / (√Σqi² × √Σri²)
    
    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        mz_tolerance: m/z tolerance for peak matching
        matched_pairs: Pre-computed matched pairs (optional)
    
    Returns:
        Cosine similarity score (0-1)
    """
    if not query_peaks or not reference_peaks:
        return 0.0
    
    if matched_pairs is None:
        matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)
    
    if not matched_pairs:
        return 0.0
    
    # Calculate dot product of matched peaks
    dot_product = sum(q.intensity * r.intensity for q, r, _ in matched_pairs)
    
    # Calculate magnitudes using all peaks
    query_magnitude = math.sqrt(sum(p.intensity ** 2 for p in query_peaks))
    ref_magnitude = math.sqrt(sum(p.intensity ** 2 for p in reference_peaks))
    
    if query_magnitude == 0 or ref_magnitude == 0:
        return 0.0
    
    return dot_product / (query_magnitude * ref_magnitude)


def calculate_dot_product(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
    matched_pairs: Optional[List[Tuple[Peak, Peak, float]]] = None,
) -> float:
    """
    Calculate normalized dot product similarity.
    
    Unlike cosine similarity, dot product considers only matched peaks
    and normalizes by the maximum possible score.
    
    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        mz_tolerance: m/z tolerance for peak matching
        matched_pairs: Pre-computed matched pairs (optional)
    
    Returns:
        Normalized dot product score (0-1)
    """
    if not query_peaks or not reference_peaks:
        return 0.0
    
    if matched_pairs is None:
        matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)
    
    if not matched_pairs:
        return 0.0
    
    dot_product = sum(q.intensity * r.intensity for q, r, _ in matched_pairs)
    
    # Normalize by maximum possible dot product
    max_possible = min(
        sum(p.intensity ** 2 for p in query_peaks),
        sum(p.intensity ** 2 for p in reference_peaks),
    )
    
    if max_possible == 0:
        return 0.0
    
    return min(1.0, dot_product / max_possible)


def calculate_weighted_cosine(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
    mz_power: float = 0.0,
    intensity_power: float = 0.5,
    matched_pairs: Optional[List[Tuple[Peak, Peak, float]]] = None,
) -> float:
    """
    Calculate weighted cosine similarity with m/z and intensity weighting.
    
    This metric applies power transformations to m/z and intensity values
    before calculating cosine similarity. Common in NIST MS Search.
    
    Formula: weight = mz^mz_power × intensity^intensity_power
    
    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        mz_tolerance: m/z tolerance for peak matching
        mz_power: Power for m/z weighting (0 = no weighting)
        intensity_power: Power for intensity weighting (0.5 = sqrt)
        matched_pairs: Pre-computed matched pairs (optional)
    
    Returns:
        Weighted cosine similarity score (0-1)
    """
    if not query_peaks or not reference_peaks:
        return 0.0
    
    if matched_pairs is None:
        matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)
    
    if not matched_pairs:
        return 0.0
    
    def weight(peak: Peak) -> float:
        """Calculate weight for a peak."""
        mz_weight = peak.mz ** mz_power if mz_power != 0 else 1.0
        int_weight = peak.intensity ** intensity_power if intensity_power != 0 else peak.intensity
        return mz_weight * int_weight
    
    # Calculate weighted dot product
    dot_product = sum(weight(q) * weight(r) for q, r, _ in matched_pairs)
    
    # Calculate weighted magnitudes
    query_magnitude = math.sqrt(sum(weight(p) ** 2 for p in query_peaks))
    ref_magnitude = math.sqrt(sum(weight(p) ** 2 for p in reference_peaks))
    
    if query_magnitude == 0 or ref_magnitude == 0:
        return 0.0
    
    return dot_product / (query_magnitude * ref_magnitude)


def calculate_spectral_entropy_similarity(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
) -> float:
    """
    Calculate spectral entropy similarity between two spectra.
    
    Spectral entropy similarity is based on information theory and measures
    how similar the information content of two spectra are. It's more robust
    to noise and small intensity variations than cosine similarity.
    
    Reference: Li et al. (2021) "Spectral entropy outperforms MS/MS dot product
    similarity for small-molecule compound identification"
    
    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        mz_tolerance: m/z tolerance for peak matching
    
    Returns:
        Spectral entropy similarity score (0-1)
    """
    if not query_peaks or not reference_peaks:
        return 0.0
    
    def normalize_intensities(peaks: List[Peak]) -> List[float]:
        """Normalize intensities to probability distribution."""
        total = sum(p.intensity for p in peaks)
        if total == 0:
            return [0.0] * len(peaks)
        return [p.intensity / total for p in peaks]
    
    def calculate_entropy(probs: List[float]) -> float:
        """Calculate Shannon entropy."""
        return -sum(p * math.log(p + 1e-10) for p in probs if p > 0)
    
    # Create combined spectrum for entropy calculation
    matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)
    
    if not matched_pairs:
        return 0.0
    
    # Build aligned intensity vectors
    query_intensities = []
    ref_intensities = []
    
    for q, r, _ in matched_pairs:
        query_intensities.append(q.intensity)
        ref_intensities.append(r.intensity)
    
    # Add unmatched peaks
    matched_query_mzs = {q.mz for q, _, _ in matched_pairs}
    matched_ref_mzs = {r.mz for _, r, _ in matched_pairs}
    
    for p in query_peaks:
        if p.mz not in matched_query_mzs:
            query_intensities.append(p.intensity)
            ref_intensities.append(0.0)
    
    for p in reference_peaks:
        if p.mz not in matched_ref_mzs:
            query_intensities.append(0.0)
            ref_intensities.append(p.intensity)
    
    # Normalize
    total_q = sum(query_intensities)
    total_r = sum(ref_intensities)
    
    if total_q == 0 or total_r == 0:
        return 0.0
    
    prob_q = [i / total_q for i in query_intensities]
    prob_r = [i / total_r for i in ref_intensities]
    
    # Calculate entropies
    entropy_q = calculate_entropy(prob_q)
    entropy_r = calculate_entropy(prob_r)
    
    # Calculate joint entropy (mixed spectrum)
    prob_mixed = [(p + q) / 2 for p, q in zip(prob_q, prob_r)]
    entropy_mixed = calculate_entropy(prob_mixed)
    
    # Spectral entropy similarity
    # 1 - (2 * H_mixed - H_q - H_r) / log(4)
    jsd = 2 * entropy_mixed - entropy_q - entropy_r
    max_jsd = 2 * math.log(2)  # Maximum Jensen-Shannon divergence
    
    similarity = 1.0 - (jsd / max_jsd)
    return max(0.0, min(1.0, similarity))


def calculate_reverse_cosine(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
) -> float:
    """
    Calculate reverse (library) cosine similarity.
    
    Reverse cosine uses only peaks present in the reference spectrum,
    ignoring additional peaks in the query. This is useful when the
    query may have additional noise or fragments.
    
    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        mz_tolerance: m/z tolerance for peak matching
    
    Returns:
        Reverse cosine similarity score (0-1)
    """
    if not query_peaks or not reference_peaks:
        return 0.0
    
    matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)
    
    if not matched_pairs:
        return 0.0
    
    # Use only matched reference peaks in calculation
    matched_ref_intensities = [r.intensity for _, r, _ in matched_pairs]
    
    # Dot product
    dot_product = sum(q.intensity * r.intensity for q, r, _ in matched_pairs)
    
    # Query magnitude (only matched peaks)
    query_magnitude = math.sqrt(sum(q.intensity ** 2 for q, _, _ in matched_pairs))
    
    # Reference magnitude (all reference peaks)
    ref_magnitude = math.sqrt(sum(p.intensity ** 2 for p in reference_peaks))
    
    if query_magnitude == 0 or ref_magnitude == 0:
        return 0.0
    
    return dot_product / (query_magnitude * ref_magnitude)


def calculate_intensity_correlation(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
) -> float:
    """
    Calculate Pearson correlation coefficient of matched peak intensities.
    
    This metric measures how well the intensity patterns of matched peaks
    correlate between the query and reference spectra.
    
    Args:
        query_peaks: Query spectrum peaks
        reference_peaks: Reference spectrum peaks
        mz_tolerance: m/z tolerance for peak matching
    
    Returns:
        Pearson correlation coefficient (-1 to 1, returned as 0-1 scale)
    """
    matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)
    
    if len(matched_pairs) < 2:
        return 0.0
    
    query_int = [q.intensity for q, _, _ in matched_pairs]
    ref_int = [r.intensity for _, r, _ in matched_pairs]
    
    # Calculate Pearson correlation
    n = len(query_int)
    mean_q = sum(query_int) / n
    mean_r = sum(ref_int) / n
    
    numerator = sum((q - mean_q) * (r - mean_r) for q, r in zip(query_int, ref_int))
    
    std_q = math.sqrt(sum((q - mean_q) ** 2 for q in query_int))
    std_r = math.sqrt(sum((r - mean_r) ** 2 for r in ref_int))
    
    if std_q == 0 or std_r == 0:
        return 0.0
    
    correlation = numerator / (std_q * std_r)
    
    # Convert from [-1, 1] to [0, 1]
    return (correlation + 1) / 2


class MetricsCalculator:
    """
    Calculator for all spectrum similarity metrics.
    
    Provides a unified interface to calculate multiple similarity metrics
    between query and reference spectra.
    """
    
    def __init__(
        self,
        mz_tolerance: float = 0.5,
        mz_power: float = 0.0,
        intensity_power: float = 0.5,
    ):
        """
        Initialize metrics calculator.
        
        Args:
            mz_tolerance: m/z tolerance for peak matching (Da)
            mz_power: Power for m/z weighting in weighted cosine
            intensity_power: Power for intensity weighting in weighted cosine
        """
        self.mz_tolerance = mz_tolerance
        self.mz_power = mz_power
        self.intensity_power = intensity_power
    
    def calculate_all_metrics(
        self,
        query_peaks: List[Peak],
        reference_peaks: List[Peak],
    ) -> SimilarityMetrics:
        """
        Calculate all similarity metrics between two spectra.
        
        Args:
            query_peaks: Query spectrum peaks
            reference_peaks: Reference spectrum peaks
        
        Returns:
            SimilarityMetrics with all calculated metrics
        """
        if not query_peaks or not reference_peaks:
            return SimilarityMetrics(
                total_query_peaks=len(query_peaks) if query_peaks else 0,
                total_reference_peaks=len(reference_peaks) if reference_peaks else 0,
            )
        
        # Match peaks once for efficiency
        matched_pairs = match_peaks(query_peaks, reference_peaks, self.mz_tolerance)
        
        # Calculate all metrics
        cosine = calculate_cosine_similarity(
            query_peaks, reference_peaks, self.mz_tolerance, matched_pairs
        )
        
        dot_product = calculate_dot_product(
            query_peaks, reference_peaks, self.mz_tolerance, matched_pairs
        )
        
        spectral_entropy = calculate_spectral_entropy_similarity(
            query_peaks, reference_peaks, self.mz_tolerance
        )
        
        weighted_cosine = calculate_weighted_cosine(
            query_peaks, reference_peaks, self.mz_tolerance,
            self.mz_power, self.intensity_power, matched_pairs
        )
        
        reverse_cosine = calculate_reverse_cosine(
            query_peaks, reference_peaks, self.mz_tolerance
        )
        
        intensity_correlation = calculate_intensity_correlation(
            query_peaks, reference_peaks, self.mz_tolerance
        )
        
        # Calculate matched peaks ratio
        matched_peaks = len(matched_pairs)
        matched_peaks_ratio = matched_peaks / max(len(query_peaks), len(reference_peaks))
        
        return SimilarityMetrics(
            cosine=cosine,
            dot_product=dot_product,
            spectral_entropy=spectral_entropy,
            weighted_cosine=weighted_cosine,
            reverse_cosine=reverse_cosine,
            matched_peaks_ratio=matched_peaks_ratio,
            intensity_correlation=intensity_correlation,
            matched_peaks=matched_peaks,
            total_query_peaks=len(query_peaks),
            total_reference_peaks=len(reference_peaks),
        )
    
    def calculate_single_metric(
        self,
        query_peaks: List[Peak],
        reference_peaks: List[Peak],
        metric: str = "cosine",
    ) -> float:
        """
        Calculate a single similarity metric.
        
        Args:
            query_peaks: Query spectrum peaks
            reference_peaks: Reference spectrum peaks
            metric: Metric name ("cosine", "dot_product", "spectral_entropy", etc.)
        
        Returns:
            Similarity score (0-1)
        """
        metric_functions = {
            "cosine": calculate_cosine_similarity,
            "dot_product": calculate_dot_product,
            "spectral_entropy": calculate_spectral_entropy_similarity,
            "weighted_cosine": lambda q, r, t: calculate_weighted_cosine(
                q, r, t, self.mz_power, self.intensity_power
            ),
            "reverse_cosine": calculate_reverse_cosine,
            "intensity_correlation": calculate_intensity_correlation,
        }
        
        if metric not in metric_functions:
            raise ValueError(f"Unknown metric: {metric}. Available: {list(metric_functions.keys())}")
        
        return metric_functions[metric](query_peaks, reference_peaks, self.mz_tolerance)


def compute_statistical_analysis(scores: List[float]) -> StatisticalAnalysis:
    """
    Compute statistical analysis of similarity scores.
    
    Args:
        scores: List of similarity scores
    
    Returns:
        StatisticalAnalysis with comprehensive statistics
    """
    if not scores:
        return StatisticalAnalysis()
    
    n = len(scores)
    scores_array = np.array(scores)
    
    mean = float(np.mean(scores_array))
    std = float(np.std(scores_array, ddof=1)) if n > 1 else 0.0
    
    # Calculate confidence interval (95%)
    if n > 1:
        se = std / math.sqrt(n)
        ci_margin = 1.96 * se
        ci_lower = mean - ci_margin
        ci_upper = mean + ci_margin
    else:
        ci_lower = ci_upper = mean
    
    return StatisticalAnalysis(
        count=n,
        mean=mean,
        std=std,
        median=float(np.median(scores_array)),
        min_value=float(np.min(scores_array)),
        max_value=float(np.max(scores_array)),
        percentile_25=float(np.percentile(scores_array, 25)),
        percentile_75=float(np.percentile(scores_array, 75)),
        percentile_90=float(np.percentile(scores_array, 90)),
        percentile_95=float(np.percentile(scores_array, 95)),
        above_threshold_50=float(np.mean(scores_array >= 0.5) * 100),
        above_threshold_70=float(np.mean(scores_array >= 0.7) * 100),
        above_threshold_90=float(np.mean(scores_array >= 0.9) * 100),
        ci_lower=max(0.0, ci_lower),
        ci_upper=min(1.0, ci_upper),
    )


def assess_match_quality(metrics: SimilarityMetrics) -> Dict[str, Any]:
    """
    Assess the quality of a spectrum match based on multiple metrics.
    
    Args:
        metrics: SimilarityMetrics from comparison
    
    Returns:
        Dictionary with quality assessment
    """
    # Quality tiers
    if metrics.cosine >= 0.9 and metrics.matched_peaks_ratio >= 0.7:
        quality_tier = "excellent"
        confidence = "high"
    elif metrics.cosine >= 0.7 and metrics.matched_peaks_ratio >= 0.5:
        quality_tier = "good"
        confidence = "medium-high"
    elif metrics.cosine >= 0.5 and metrics.matched_peaks_ratio >= 0.3:
        quality_tier = "moderate"
        confidence = "medium"
    elif metrics.cosine >= 0.3:
        quality_tier = "poor"
        confidence = "low"
    else:
        quality_tier = "no_match"
        confidence = "none"
    
    # Check for metric consistency
    metric_values = [
        metrics.cosine,
        metrics.dot_product,
        metrics.spectral_entropy,
    ]
    metric_std = float(np.std(metric_values)) if len(metric_values) > 1 else 0.0
    metrics_consistent = metric_std < 0.15
    
    return {
        "quality_tier": quality_tier,
        "confidence": confidence,
        "metrics_consistent": metrics_consistent,
        "metric_std": metric_std,
        "primary_metric": metrics.cosine,
        "recommendation": _get_quality_recommendation(quality_tier, metrics),
    }


def _get_quality_recommendation(quality_tier: str, metrics: SimilarityMetrics) -> str:
    """Generate recommendation based on quality tier."""
    recommendations = {
        "excellent": "Strong match. Consider validated identification.",
        "good": "Good match. Additional verification recommended.",
        "moderate": "Moderate match. Manual review recommended.",
        "poor": "Poor match. Consider alternative candidates.",
        "no_match": "No significant match found.",
    }
    return recommendations.get(quality_tier, "Unknown quality tier.")


# ============================================================================
# Specialized Metrics for QCxMS2 Validation
# ============================================================================

def calculate_fragmentation_pattern_match(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
    top_n: int = 10,
) -> float:
    """
    Calculate fragmentation pattern match focusing on top peaks.
    
    This metric is designed for comparing calculated (QCxMS2) spectra with
    experimental data. It focuses on whether the calculation correctly
    predicts the PRESENCE of major fragments, with less emphasis on exact
    intensity ratios (which vary with experimental conditions).
    
    Algorithm:
    1. Identify top-N peaks by intensity in both spectra
    2. Count how many of reference's top-N peaks are found in query
    3. Weight matches by reference intensity rank
    
    Args:
        query_peaks: Query spectrum peaks (calculated)
        reference_peaks: Reference spectrum peaks (experimental)
        mz_tolerance: m/z tolerance for peak matching (Da)
        top_n: Number of top peaks to consider
    
    Returns:
        Pattern match score (0-1)
    
    Notes:
        - Score of 1.0 means all top-N reference peaks were found in query
        - This metric is asymmetric: it asks "does the calculation predict
          the observed fragments?" not "are all calculated fragments observed?"
    """
    if not query_peaks or not reference_peaks:
        return 0.0
    
    # Sort reference peaks by intensity (descending)
    ref_sorted = sorted(reference_peaks, key=lambda p: p.intensity, reverse=True)
    top_ref = ref_sorted[:top_n]
    
    # Check how many top reference peaks are found in query
    matched = 0
    total_weight = 0
    
    for rank, ref_peak in enumerate(top_ref):
        # Weight by rank (first peak = top_n points, last = 1 point)
        weight = top_n - rank
        total_weight += weight
        
        # Check if this reference peak is found in query
        for q_peak in query_peaks:
            if abs(q_peak.mz - ref_peak.mz) <= mz_tolerance:
                matched += weight
                break
    
    return matched / total_weight if total_weight > 0 else 0.0


def calculate_mass_accuracy_score(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
) -> Tuple[float, List[float]]:
    """
    Calculate mass accuracy score and report m/z errors.
    
    This metric specifically assesses how accurately QCxMS2 predicts
    m/z values. Useful for validating computational methodology.
    
    Args:
        query_peaks: Query spectrum peaks (calculated)
        reference_peaks: Reference spectrum peaks (experimental)
        mz_tolerance: m/z tolerance for peak matching (Da)
    
    Returns:
        Tuple of (accuracy_score, list_of_mz_errors)
        - accuracy_score: Mean inverse error (higher = better accuracy)
        - list_of_mz_errors: Actual m/z differences for matched peaks
    """
    if not query_peaks or not reference_peaks:
        return 0.0, []
    
    matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)
    
    if not matched_pairs:
        return 0.0, []
    
    mz_errors = [abs(q.mz - r.mz) for q, r, _ in matched_pairs]
    
    # Calculate accuracy score as mean inverse error
    # Use 0.001 as floor to avoid division by zero
    mean_error = sum(mz_errors) / len(mz_errors)
    accuracy_score = 1.0 / (1.0 + mean_error * 100)  # Scale: 0.01 Da error -> 0.5 score
    
    return accuracy_score, mz_errors


def calculate_intensity_rank_correlation(
    query_peaks: List[Peak],
    reference_peaks: List[Peak],
    mz_tolerance: float = 0.5,
) -> float:
    """
    Calculate Spearman rank correlation of matched peak intensities.
    
    This metric checks if the RELATIVE ordering of fragment intensities
    matches between calculated and experimental spectra, even if the
    absolute intensity ratios differ.
    
    Args:
        query_peaks: Query spectrum peaks (calculated)
        reference_peaks: Reference spectrum peaks (experimental)
        mz_tolerance: m/z tolerance for peak matching (Da)
    
    Returns:
        Rank correlation score (-1 to 1, scaled to 0-1)
    
    Notes:
        - Score of 1.0 means perfect rank agreement
        - Score of 0.5 means no correlation
        - This is useful when intensity ratios differ but relative
          ordering should be preserved
    """
    if not query_peaks or not reference_peaks:
        return 0.0
    
    matched_pairs = match_peaks(query_peaks, reference_peaks, mz_tolerance)
    
    if len(matched_pairs) < 3:
        return 0.0  # Need at least 3 points for meaningful correlation
    
    # Get intensity ranks
    query_intensities = [q.intensity for q, _, _ in matched_pairs]
    ref_intensities = [r.intensity for _, r, _ in matched_pairs]
    
    # Calculate ranks
    def get_ranks(values):
        sorted_idx = sorted(range(len(values)), key=lambda i: values[i])
        ranks = [0] * len(values)
        for rank, idx in enumerate(sorted_idx):
            ranks[idx] = rank + 1
        return ranks
    
    q_ranks = get_ranks(query_intensities)
    r_ranks = get_ranks(ref_intensities)
    
    # Spearman rank correlation
    n = len(q_ranks)
    d_squared_sum = sum((q - r) ** 2 for q, r in zip(q_ranks, r_ranks))
    rho = 1 - (6 * d_squared_sum) / (n * (n ** 2 - 1))
    
    # Scale from [-1, 1] to [0, 1]
    return (rho + 1) / 2


@dataclass
class QCxMS2ValidationMetrics:
    """
    Specialized metrics for validating QCxMS2 calculations against experiment.
    
    These metrics are designed to address the known characteristics of
    computational mass spectra:
    - Excellent m/z accuracy (usually < 0.01 Da)
    - Variable intensity ratios (due to sampling statistics)
    - Possible prediction of minor peaks not observed experimentally
    """
    
    # Standard metrics
    cosine_similarity: float = 0.0
    weighted_cosine: float = 0.0
    
    # QCxMS2-specific metrics
    fragmentation_pattern_match: float = 0.0  # Top-N peak detection
    mass_accuracy_score: float = 0.0          # m/z precision
    intensity_rank_correlation: float = 0.0   # Relative ordering
    
    # Peak matching stats
    matched_peaks: int = 0
    query_peaks: int = 0
    reference_peaks: int = 0
    mean_mz_error: float = 0.0
    max_mz_error: float = 0.0
    
    # Derived scores
    @property
    def overall_validation_score(self) -> float:
        """
        Combined validation score emphasizing fragmentation pattern match.
        
        Weights:
        - 40% Fragmentation pattern match (did we find the right fragments?)
        - 30% Mass accuracy (are the m/z values correct?)
        - 20% Intensity rank correlation (is the relative ordering right?)
        - 10% Cosine similarity (traditional metric)
        """
        return (
            0.40 * self.fragmentation_pattern_match +
            0.30 * self.mass_accuracy_score +
            0.20 * self.intensity_rank_correlation +
            0.10 * self.cosine_similarity
        )
    
    @property
    def quality_assessment(self) -> str:
        """Assess validation quality."""
        score = self.overall_validation_score
        if score >= 0.85:
            return "Excellent - QCxMS2 accurately predicts fragmentation"
        elif score >= 0.70:
            return "Good - Major fragmentation pathways captured"
        elif score >= 0.50:
            return "Moderate - Some fragmentation pathways captured"
        else:
            return "Poor - Significant disagreement with experiment"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "cosine_similarity": self.cosine_similarity,
            "weighted_cosine": self.weighted_cosine,
            "fragmentation_pattern_match": self.fragmentation_pattern_match,
            "mass_accuracy_score": self.mass_accuracy_score,
            "intensity_rank_correlation": self.intensity_rank_correlation,
            "overall_validation_score": self.overall_validation_score,
            "quality_assessment": self.quality_assessment,
            "matched_peaks": self.matched_peaks,
            "query_peaks": self.query_peaks,
            "reference_peaks": self.reference_peaks,
            "mean_mz_error": self.mean_mz_error,
            "max_mz_error": self.max_mz_error,
        }


def calculate_qcxms2_validation_metrics(
    calculated_peaks: List[Peak],
    experimental_peaks: List[Peak],
    mz_tolerance: float = 0.5,
    top_n_peaks: int = 10,
) -> QCxMS2ValidationMetrics:
    """
    Calculate all QCxMS2 validation metrics.
    
    This is the main function for validating QCxMS2 calculations against
    experimental mass spectra.
    
    Args:
        calculated_peaks: Peaks from QCxMS2 calculation
        experimental_peaks: Peaks from experimental spectrum (MassBank, etc.)
        mz_tolerance: m/z tolerance for peak matching (Da)
        top_n_peaks: Number of top peaks for pattern matching
    
    Returns:
        QCxMS2ValidationMetrics with all metrics calculated
    
    Example:
        >>> from database_integration.benchmark.metrics import (
        ...     calculate_qcxms2_validation_metrics, Peak
        ... )
        >>> calc = [Peak(78.05, 100), Peak(39.02, 40), Peak(51.02, 15)]
        >>> expt = [Peak(78.05, 100), Peak(51.02, 20), Peak(39.02, 10)]
        >>> metrics = calculate_qcxms2_validation_metrics(calc, expt)
        >>> print(f"Validation: {metrics.overall_validation_score:.1%}")
    """
    if not calculated_peaks or not experimental_peaks:
        return QCxMS2ValidationMetrics(
            query_peaks=len(calculated_peaks) if calculated_peaks else 0,
            reference_peaks=len(experimental_peaks) if experimental_peaks else 0,
        )
    
    # Standard metrics
    matched_pairs = match_peaks(calculated_peaks, experimental_peaks, mz_tolerance)
    cosine = calculate_cosine_similarity(
        calculated_peaks, experimental_peaks, mz_tolerance, matched_pairs
    )
    weighted = calculate_weighted_cosine(
        calculated_peaks, experimental_peaks, mz_tolerance, 
        mz_power=0.0, intensity_power=0.5, matched_pairs=matched_pairs
    )
    
    # QCxMS2-specific metrics
    pattern_match = calculate_fragmentation_pattern_match(
        calculated_peaks, experimental_peaks, mz_tolerance, top_n_peaks
    )
    
    accuracy, mz_errors = calculate_mass_accuracy_score(
        calculated_peaks, experimental_peaks, mz_tolerance
    )
    
    rank_correlation = calculate_intensity_rank_correlation(
        calculated_peaks, experimental_peaks, mz_tolerance
    )
    
    return QCxMS2ValidationMetrics(
        cosine_similarity=cosine,
        weighted_cosine=weighted,
        fragmentation_pattern_match=pattern_match,
        mass_accuracy_score=accuracy,
        intensity_rank_correlation=rank_correlation,
        matched_peaks=len(matched_pairs),
        query_peaks=len(calculated_peaks),
        reference_peaks=len(experimental_peaks),
        mean_mz_error=sum(mz_errors) / len(mz_errors) if mz_errors else 0.0,
        max_mz_error=max(mz_errors) if mz_errors else 0.0,
    )
