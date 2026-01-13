#!/usr/bin/env python3
"""
Example: Parse QCxMS2 output and compare with experimental spectra.

This example demonstrates how to:
1. Parse QCxMS2 calculation output files
2. Convert to Peak objects for database integration
3. Compare calculated spectrum with experimental data from MassBank

Usage:
    python parse_qcxms2_output.py /path/to/qcxms2/calculation/directory

Requirements:
    - A completed QCxMS2 calculation directory with peaks.dat
    - Optional: Network access for MassBank comparison
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from database_integration import (
    # QCxMS2 parsers
    parse_qcxms2_results,
    parse_qcxms2_output_dir,
    QCxMS2Result,
    # Database connectors
    create_ms_connector,
    DatabaseType,
    SpectrumSearchCriteria,
    # Comparison utilities  
    compare_spectra,
    Peak,
)


def example_parse_peaks_file():
    """Example: Parse a single peaks.dat file."""
    print("=" * 60)
    print("Example 1: Parse peaks.dat file")
    print("=" * 60)
    
    # Create sample spectrum text (simulating peaks.dat content)
    sample_spectrum = """
# QCxMS2 calculated spectrum for ethanol
# m/z  intensity (normalized to 10000)
31.018538  10000.0
46.041865  1500.0
45.033290  500.0
29.002740  3000.0
28.031299  1000.0
27.023425  200.0
"""
    
    from database_integration.parsers import parse_qcxms2_spectrum
    
    # Parse spectrum text
    peaks = parse_qcxms2_spectrum(sample_spectrum, normalize_to=100.0)
    
    print(f"\nParsed {len(peaks)} peaks:")
    print("-" * 40)
    print(f"{'m/z':>12}  {'Intensity':>12}")
    print("-" * 40)
    for peak in peaks:
        print(f"{peak.mz:>12.4f}  {peak.intensity:>12.2f}")
    
    # Find base peak
    base_peak = max(peaks, key=lambda p: p.intensity)
    print(f"\nBase peak: m/z = {base_peak.mz:.4f} (100%)")
    
    return peaks


def example_parse_calculation_dir(calc_dir: str = None):
    """Example: Parse a complete QCxMS2 calculation directory."""
    print("\n" + "=" * 60)
    print("Example 2: Parse complete QCxMS2 calculation directory")
    print("=" * 60)
    
    if calc_dir is None:
        print("\nNo calculation directory provided.")
        print("Usage: python parse_qcxms2_output.py /path/to/qcxms2/calc/")
        print("\nExpected files in directory:")
        print("  - peaks.dat    : Main spectrum (m/z, intensity)")
        print("  - allpeaks.dat : Fragment assignments")
        print("  - qcxms2.out   : Calculation output with metadata")
        return None
    
    calc_path = Path(calc_dir)
    if not calc_path.exists():
        print(f"\nError: Directory not found: {calc_path}")
        return None
    
    try:
        result = parse_qcxms2_output_dir(calc_path, normalize_to=100.0)
        
        print(f"\nCalculation directory: {result.source_dir}")
        print(f"Ionization mode: {result.mode.value}")
        print(f"QCxMS2 version: {result.version}")
        print(f"Input file: {result.input_file}")
        print(f"Calculation level: {result.calculation_level}")
        if result.ip_ev:
            print(f"Ionization potential: {result.ip_ev:.3f} eV")
        
        print(f"\nSpectrum: {len(result.peaks)} peaks")
        print("-" * 50)
        print(f"{'m/z':>12}  {'Intensity':>12}  {'Annotation'}")
        print("-" * 50)
        for peak in result.peaks[:15]:  # Show first 15 peaks
            ann = peak.annotation or ""
            print(f"{peak.mz:>12.4f}  {peak.intensity:>12.2f}  {ann}")
        if len(result.peaks) > 15:
            print(f"  ... and {len(result.peaks) - 15} more peaks")
        
        base_peak = result.get_base_peak()
        if base_peak:
            print(f"\nBase peak: m/z = {base_peak.mz:.4f}")
        
        mol_ion = result.get_molecular_ion_peak()
        if mol_ion:
            print(f"Molecular ion (M+): m/z = {mol_ion.mz:.4f} ({mol_ion.intensity:.1f}%)")
        
        if result.fragments:
            print(f"\nFragment assignments: {len(result.fragments)} entries")
            print("-" * 60)
            print(f"{'Path':<20} {'m/z':>10} {'I (%)':>8} {'Formula'}")
            print("-" * 60)
            for frag in result.fragments[:10]:
                formula = frag.sum_formula or ""
                print(f"{frag.fragment_path:<20} {frag.mz:>10.4f} {frag.intensity:>8.2f} {formula}")
        
        return result
        
    except Exception as e:
        print(f"\nError parsing calculation: {e}")
        return None


def example_compare_with_massbank(qcxms2_peaks: list = None):
    """Example: Compare QCxMS2 spectrum with MassBank experimental data."""
    print("\n" + "=" * 60)
    print("Example 3: Compare with MassBank experimental spectrum")
    print("=" * 60)
    
    # Use provided peaks or create sample ethanol spectrum
    if qcxms2_peaks is None:
        qcxms2_peaks = [
            Peak(mz=31.0, intensity=100.0, annotation="CH3O+"),
            Peak(mz=46.0, intensity=15.0, annotation="M+"),
            Peak(mz=45.0, intensity=5.0, annotation="C2H5O+"),
            Peak(mz=29.0, intensity=29.0, annotation="CHO+"),
            Peak(mz=28.0, intensity=11.0, annotation="C2H4+"),
            Peak(mz=27.0, intensity=2.0, annotation="C2H3+"),
        ]
        print("\nUsing sample ethanol spectrum (QCxMS2-like)")
    
    print(f"\nQCxMS2 calculated spectrum: {len(qcxms2_peaks)} peaks")
    
    # Try to connect to MassBank
    try:
        print("\nConnecting to MassBank EU...")
        massbank = create_ms_connector(DatabaseType.MASSBANK)
        
        # Search for ethanol EI spectra
        criteria = SpectrumSearchCriteria(
            formula="C2H6O",
            ionization_mode="EI",
            max_results=5,
        )
        
        result = massbank.search_spectrum(criteria)
        
        if result.success and result.spectra:
            print(f"Found {len(result.spectra)} experimental spectra")
            
            # Compare with first experimental spectrum
            exp_spectrum = result.spectra[0]
            print(f"\nComparing with: {exp_spectrum.name}")
            print(f"Source: {exp_spectrum.source}")
            
            comparison = compare_spectra(
                query_peaks=qcxms2_peaks,
                reference_peaks=exp_spectrum.peaks,
                mz_tolerance=0.5,
                query_id="QCxMS2_calculated",
                reference_id=exp_spectrum.spectrum_id,
            )
            
            print(f"\nComparison results:")
            print(f"  Similarity score: {comparison.similarity_score:.3f}")
            print(f"  Matched peaks: {comparison.matched_peaks} / {comparison.total_peaks_query}")
            print(f"  m/z tolerance: {comparison.mz_tolerance} Da")
            
            if comparison.matched_peak_pairs:
                print(f"\nMatched peak pairs:")
                print(f"  {'Calc m/z':>10} {'Exp m/z':>10} {'Calc I':>8} {'Exp I':>8}")
                for calc_peak, exp_peak in comparison.matched_peak_pairs[:10]:
                    print(f"  {calc_peak.mz:>10.2f} {exp_peak.mz:>10.2f} "
                          f"{calc_peak.intensity:>8.1f} {exp_peak.intensity:>8.1f}")
        else:
            print(f"No experimental spectra found: {result.error_message}")
            
    except Exception as e:
        print(f"\nMassBank comparison failed: {e}")
        print("(This is expected if running without network access)")


def example_export_to_common_formats(peaks: list):
    """Example: Export spectrum to common formats."""
    print("\n" + "=" * 60)
    print("Example 4: Export to common formats")
    print("=" * 60)
    
    if not peaks:
        print("No peaks to export")
        return
    
    # MSP format (NIST)
    print("\nMSP format (NIST):")
    print("-" * 40)
    print("Name: QCxMS2 Calculated Spectrum")
    print("Num Peaks: {}".format(len(peaks)))
    for peak in peaks[:5]:
        print(f"{peak.mz:.4f} {peak.intensity:.1f}")
    if len(peaks) > 5:
        print(f"... ({len(peaks) - 5} more peaks)")
    
    # CSV format
    print("\n\nCSV format:")
    print("-" * 40)
    print("m/z,intensity,annotation")
    for peak in peaks[:5]:
        ann = peak.annotation or ""
        print(f"{peak.mz:.4f},{peak.intensity:.2f},{ann}")
    if len(peaks) > 5:
        print(f"... ({len(peaks) - 5} more peaks)")


def main():
    """Run all examples."""
    print("QCxMS2 Output Parser Examples")
    print("=" * 60)
    
    # Example 1: Parse spectrum text
    peaks = example_parse_peaks_file()
    
    # Example 2: Parse calculation directory (if provided)
    calc_dir = sys.argv[1] if len(sys.argv) > 1 else None
    result = example_parse_calculation_dir(calc_dir)
    
    # Use parsed peaks if available, otherwise use sample
    if result and result.peaks:
        peaks_for_comparison = result.peaks
    else:
        peaks_for_comparison = peaks
    
    # Example 3: Compare with MassBank
    example_compare_with_massbank(peaks_for_comparison)
    
    # Example 4: Export formats
    example_export_to_common_formats(peaks_for_comparison)
    
    print("\n" + "=" * 60)
    print("Examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
