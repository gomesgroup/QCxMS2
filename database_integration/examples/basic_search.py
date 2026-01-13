#!/usr/bin/env python3
"""
Basic example: Search MassBank for glucose EI spectra.

Usage:
    python basic_search.py
"""

import sys
from pathlib import Path

# Add database_integration to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from database_integration import (
    DatabaseType,
    IonizationMode,
    SpectrumSearchCriteria,
    create_ms_connector,
)


def main():
    print("=== QCxMS2 Database Integration: Basic Search Example ===\n")

    # Create MassBank connector
    print("Connecting to MassBank EU...")
    massbank = create_ms_connector(DatabaseType.MASSBANK)

    # Test connection
    if not massbank.test_connection():
        print("❌ Failed to connect to MassBank")
        return

    print("✅ Connected to MassBank\n")

    # Search for glucose EI spectra
    print("Searching for C6H12O6 (glucose) with EI ionization...")
    criteria = SpectrumSearchCriteria(
        formula="C6H12O6",
        ionization_mode=IonizationMode.EI,
        max_results=10,
    )

    results = massbank.search_spectrum(criteria)

    if not results.success:
        print(f"❌ Search failed: {results.error_message}")
        return

    print(f"✅ Found {results.total_results} spectra\n")

    # Display results
    for i, spectrum in enumerate(results.spectra, 1):
        print(f"{i}. {spectrum.spectrum_id}")
        print(f"   Name: {spectrum.name}")
        print(f"   Formula: {spectrum.formula}")
        print(f"   Exact Mass: {spectrum.exact_mass:.4f}" if spectrum.exact_mass else "")
        print(f"   Peaks: {len(spectrum.peaks)}")
        print(f"   Instrument: {spectrum.instrument}")
        print(f"   URL: {spectrum.url}")
        print()

    # Get detailed info for first spectrum
    if results.spectra:
        first_spectrum = results.spectra[0]
        print(f"\nDetailed info for {first_spectrum.spectrum_id}:")
        print(f"  InChI Key: {first_spectrum.inchi_key}")
        print(f"  Collision Energy: {first_spectrum.collision_energy}")

        print(f"\n  First 10 peaks:")
        for peak in first_spectrum.peaks[:10]:
            print(f"    m/z {peak.mz:.2f}: {peak.intensity:.1f}")

    print(f"\nCitation: {results.citation}")


if __name__ == "__main__":
    main()
