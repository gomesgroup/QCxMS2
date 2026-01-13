#!/usr/bin/env python3
"""
Example: Searching MassSpecGym benchmark dataset

MassSpecGym is a NeurIPS 2024 Spotlight benchmark containing 231K curated
MS/MS spectra from 29K molecules, combining MassBank, MoNA, and GNPS data.

This example shows how to:
1. Load the MassSpecGym dataset
2. Search by molecular formula
3. Search by InChI key
4. Get dataset statistics
5. Iterate over spectra
"""

import sys
sys.path.insert(0, "/mnt/beegfs/software/qcxms2")

from database_integration import (
    create_ms_connector,
    DatabaseType,
    SpectrumSearchCriteria,
    IonizationMode,
)
from database_integration.factory import create_massspecgym_connector


def main():
    """Demonstrate MassSpecGym database integration."""

    print("=" * 70)
    print("MassSpecGym Database Integration Example")
    print("NeurIPS 2024 Spotlight - 231K curated MS/MS spectra")
    print("=" * 70)
    print()

    # Method 1: Using factory function directly
    print("Creating MassSpecGym connector...")
    print("(First run downloads ~2GB dataset from Hugging Face)")
    print()

    connector = create_massspecgym_connector()

    # Test connection (loads dataset)
    print("Testing connection / loading dataset...")
    if not connector.test_connection():
        print("ERROR: Failed to load MassSpecGym dataset")
        print("Make sure 'datasets' package is installed: pip install datasets")
        return 1

    # Get dataset statistics
    print("\n--- Dataset Statistics ---")
    stats = connector.get_dataset_statistics()
    for key, value in stats.items():
        print(f"  {key}: {value}")
    print()

    # Example 1: Search by molecular formula
    print("--- Example 1: Search by Molecular Formula ---")
    criteria = SpectrumSearchCriteria(
        formula="C8H10N4O2",  # Caffeine
        max_results=5,
    )
    results = connector.search_spectrum(criteria)

    if results.success:
        print(f"Found {results.total_results} spectra for C8H10N4O2 (caffeine)")
        for i, spectrum in enumerate(results.spectra[:3]):
            print(f"\n  [{i+1}] {spectrum.spectrum_id}")
            print(f"      SMILES: {spectrum.smiles[:50] if spectrum.smiles else 'N/A'}...")
            print(f"      InChI Key: {spectrum.inchi_key}")
            print(f"      Precursor m/z: {spectrum.precursor_mz}")
            print(f"      Adduct: {spectrum.precursor_type}")
            print(f"      Peaks: {len(spectrum.peaks)}")
    else:
        print(f"Search failed: {results.error_message}")
    print()

    # Example 2: Search by InChI Key
    print("--- Example 2: Search by InChI Key ---")
    # Aspirin InChI Key (2D portion)
    criteria = SpectrumSearchCriteria(
        inchi_key="BSYNRYMUTXBXSQ",  # Aspirin 2D InChI key prefix
        max_results=5,
    )
    results = connector.search_spectrum(criteria)

    if results.success:
        print(f"Found {results.total_results} spectra for InChI key BSYNRYMUTXBXSQ...")
        for spectrum in results.spectra[:2]:
            print(f"\n  ID: {spectrum.spectrum_id}")
            print(f"  Formula: {spectrum.formula}")
            print(f"  Peaks: {len(spectrum.peaks)}")
    else:
        print(f"Search failed: {results.error_message}")
    print()

    # Example 3: Get spectrum by index
    print("--- Example 3: Get Spectrum by Index ---")
    spectrum = connector.get_spectrum_by_id("0")  # First spectrum
    if spectrum:
        print(f"Spectrum 0:")
        print(f"  ID: {spectrum.spectrum_id}")
        print(f"  SMILES: {spectrum.smiles[:60] if spectrum.smiles else 'N/A'}...")
        print(f"  Formula: {spectrum.formula}")
        print(f"  Precursor m/z: {spectrum.precursor_mz}")
        print(f"  Number of peaks: {len(spectrum.peaks)}")
        if spectrum.peaks:
            print(f"  Top 5 peaks by intensity:")
            sorted_peaks = sorted(spectrum.peaks, key=lambda p: p.intensity, reverse=True)[:5]
            for p in sorted_peaks:
                print(f"    m/z {p.mz:.2f}: {p.intensity:.1f}")
    print()

    # Example 4: Iterate over spectra (first 5)
    print("--- Example 4: Iterate Over Spectra ---")
    print("Showing first 5 spectra:")
    for i, spectrum in enumerate(connector.iterate_spectra(max_spectra=5)):
        print(f"  [{i}] {spectrum.formula or 'N/A':15} | {len(spectrum.peaks):4} peaks | {spectrum.smiles[:30] if spectrum.smiles else 'N/A'}...")
    print()

    # Print citation
    print("--- Citation ---")
    print(connector.get_citation())
    print()

    print("=" * 70)
    print("Example completed successfully!")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
