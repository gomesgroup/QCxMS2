"""
QCxMS2 output file parsers.

This module provides parsers for reading QCxMS2 calculation output files
and converting them to Peak objects for comparison with experimental spectra.
"""

from .qcxms2_parser import (
    QCxMS2Result,
    QCxMS2FragmentInfo,
    parse_qcxms2_results,
    parse_qcxms2_spectrum,
    parse_qcxms2_allpeaks,
    parse_qcxms2_output_dir,
)

__all__ = [
    "QCxMS2Result",
    "QCxMS2FragmentInfo",
    "parse_qcxms2_results",
    "parse_qcxms2_spectrum",
    "parse_qcxms2_allpeaks",
    "parse_qcxms2_output_dir",
]
