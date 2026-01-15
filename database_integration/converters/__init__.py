"""
QCxMS2 Format Converters

Convert between QCxMS2 output formats and standard mass spectrometry formats.

Key Converters:
- peaks_to_mgf: Convert QCxMS2 peaks.dat to MGF format (for DreaMS)
- peaks_to_msp: Convert QCxMS2 peaks.dat to MSP format (NIST format)
- mgf_to_peaks: Convert MGF back to peaks format
"""

from .peaks_to_mgf import (
    peaks_to_mgf,
    peaks_to_mgf_batch,
    PeaksToMGFConverter,
)

__all__ = [
    "peaks_to_mgf",
    "peaks_to_mgf_batch",
    "PeaksToMGFConverter",
]
