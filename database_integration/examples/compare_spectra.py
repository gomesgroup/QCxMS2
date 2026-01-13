#!/usr/bin/env python3
"""
Example: Compare calculated QCxMS2 spectrum with experimental database.

Usage:
    python compare_spectra.py
"""

import sys
from pathlib import Path

# Add database_integration to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from database_integration import (
    DatabaseType,
    IonizationMode,
    Peak,
    SimilarityMetric,
    SpectrumSearchCriteria,
    compare_spectra,
    create_ms_connector,
)


def main():
    print("=== QCxMS2 Database Integration: Spectrum Comparison Example ===\n")

    # Simulated QCxMS2 calculated peaks for acetone (C3H6O)
    print("Calculated QCxMS2 spectrum for acetone (C3H6O):")
    calculated_peaks = [
        Peak(mz=15.0, intensity=15.0),
        Peak(mz=26.0, intensity=10.0),
        Peak(mz=27.0, intensity=25.0),
        Peak(mz=28.0, intensity=20.0),
        Peak(mz=29.0, intensity=15.0),
        Peak(mz=42.0, intensity=50.0),
        Peak(mz=43.0, intensity=100.0),   # Base peak
        Peak(mz=58.0, intensity=85.0),     # M+
    ]

    for peak in calculated_peaks:
        print(f"  m/z {peak.mz:.1f}: {peak.intensity:.1f}")

    # Search MassBank for acetone EI spectra
    print("\n" + "=" * 60)
    print("Searching MassBank for experimental acetone EI spectra...\n")

    massbank = create_ms_connector(DatabaseType.MASSBANK)

    criteria = SpectrumSearchCriteria(
        formula="C3H6O",
        ionization_mode=IonizationMode.EI,
        max_results=20,
    )

    results = massbank.search_spectrum(criteria)

    if not results.success:
        print(f"❌ Search failed: {results.error_message}")
        return

    print(f"✅ Found {results.total_results} experimental spectra\n")

    # Compare calculated spectrum with each experimental spectrum
    print("=" * 60)
    print("Comparing calculated vs. experimental spectra:\n")

    comparisons = []
    for exp_spectrum in results.spectra:
        comparison = compare_spectra(
            calculated_peaks,
            exp_spectrum.peaks,
            metric=SimilarityMetric.COSINE,
            mz_tolerance=0.5,
            query_id="QCxMS2_calculated",
            reference_id=exp_spectrum.spectrum_id,
        )

        comparisons.append((exp_spectrum, comparison))

    # Sort by similarity score
    comparisons.sort(key=lambda x: x[1].similarity_score, reverse=True)

    # Display top 10 matches
    print("Top 10 matches (by cosine similarity):\n")
    for i, (exp_spectrum, comparison) in enumerate(comparisons[:10], 1):
        print(f"{i}. {exp_spectrum.spectrum_id} ({exp_spectrum.name})")
        print(f"   Similarity: {comparison.similarity_score:.4f}")
        print(f"   Matched peaks: {comparison.matched_peaks}/{comparison.total_peaks_query}")
        print(f"   Instrument: {exp_spectrum.instrument}")
        print(f"   URL: {exp_spectrum.url}")
        print()

    # Detailed comparison with best match
    if comparisons:
        best_spectrum, best_comparison = comparisons[0]
        print("=" * 60)
        print(f"Detailed comparison with best match: {best_spectrum.name}")
        print(f"Similarity score: {best_comparison.similarity_score:.4f}\n")

        print("Matched peaks:")
        for calc_peak, exp_peak in best_comparison.matched_peak_pairs[:15]:
            mz_diff = abs(calc_peak.mz - exp_peak.mz)
            int_diff = abs(calc_peak.intensity - exp_peak.intensity)
            print(
                f"  Calc m/z {calc_peak.mz:.2f} ({calc_peak.intensity:.1f}) "
                f"↔ Exp m/z {exp_peak.mz:.2f} ({exp_peak.intensity:.1f}) "
                f"[Δm/z={mz_diff:.3f}, ΔI={int_diff:.1f}]"
            )

        # Try other similarity metrics
        print("\n" + "=" * 60)
        print("Comparison using different similarity metrics:\n")

        for metric in [SimilarityMetric.COSINE, SimilarityMetric.DOT_PRODUCT, SimilarityMetric.EUCLIDEAN]:
            comp = compare_spectra(
                calculated_peaks,
                best_spectrum.peaks,
                metric=metric,
                mz_tolerance=0.5,
            )
            print(f"  {metric.value:20s}: {comp.similarity_score:.4f}")

    print("\n" + "=" * 60)
    print("✅ Comparison complete!\n")


if __name__ == "__main__":
    main()
