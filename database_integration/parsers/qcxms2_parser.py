"""
QCxMS2 output file parser.

This module provides functions to parse QCxMS2 mass spectrometry calculation outputs
and convert them to Peak objects for comparison with experimental spectra.

QCxMS2 produces several output files:
- peaks.dat: Main spectrum file (m/z and intensity, space-separated)
- peaks.csv: Same data in CSV format
- allpeaks.dat: Fragment assignments with m/z and intensities
- fragments file: Contains fragment-specific information

Output intensity normalization:
- QCxMS2 normalizes intensities to 10000 (NIST convention)
- This parser normalizes to 100 by default for compatibility with database_integration

References:
    - J. Gorges, S. Grimme, PCCP 2025, 27, 6899-6911 (QCxMS2 method)
    - J. Gorges, M. Engeser, S. Grimme, JASMS 2025 (CID mode)
"""

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from ..data_classes import Peak
from ..enums import IonizationMode, SpectrumType

logger = logging.getLogger(__name__)


@dataclass
class QCxMS2FragmentInfo:
    """Information about a fragment in the QCxMS2 calculation.
    
    Attributes:
        fragment_path: Path to fragment directory (e.g., "p1f1", "p2f2")
        mz: Mass-to-charge ratio
        intensity: Relative intensity (normalized)
        sum_formula: Molecular formula of the fragment (e.g., "C2H5O")
        fragment_type: Type of fragment ("fragmentpair", "isomer", "input structure")
        reaction_energy: Reaction energy in kcal/mol (if available)
        barrier: Barrier height in kcal/mol (if available)
    """
    fragment_path: str
    mz: float
    intensity: float
    sum_formula: Optional[str] = None
    fragment_type: Optional[str] = None
    reaction_energy: Optional[float] = None
    barrier: Optional[float] = None
    
    def to_peak(self, annotation: Optional[str] = None) -> Peak:
        """Convert to Peak object.
        
        Args:
            annotation: Optional custom annotation. If None, uses sum_formula.
            
        Returns:
            Peak object with mz, intensity, and annotation.
        """
        if annotation is None:
            annotation = self.sum_formula
        return Peak(mz=self.mz, intensity=self.intensity, annotation=annotation)


@dataclass
class QCxMS2Result:
    """Complete QCxMS2 calculation result.
    
    This class holds all information from a QCxMS2 calculation including
    the spectrum peaks, fragment information, and calculation metadata.
    
    Attributes:
        peaks: List of Peak objects representing the mass spectrum
        fragments: List of QCxMS2FragmentInfo with detailed fragment data
        mode: Ionization mode (EI, CID, or DEA)
        input_file: Original input file name
        molecular_formula: Molecular formula of the input molecule
        ip_ev: Ionization potential in eV
        calculation_level: QM method used (e.g., "gfn2")
        version: QCxMS2 version string
        source_dir: Directory containing the calculation
        metadata: Additional metadata from the calculation
    """
    peaks: List[Peak] = field(default_factory=list)
    fragments: List[QCxMS2FragmentInfo] = field(default_factory=list)
    mode: IonizationMode = IonizationMode.EI
    input_file: Optional[str] = None
    molecular_formula: Optional[str] = None
    ip_ev: Optional[float] = None
    calculation_level: Optional[str] = None
    version: Optional[str] = None
    source_dir: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Sort peaks by m/z after initialization."""
        if self.peaks:
            self.peaks = sorted(self.peaks, key=lambda p: p.mz)
    
    def get_base_peak(self) -> Optional[Peak]:
        """Get the base peak (highest intensity).
        
        Returns:
            Peak with highest intensity, or None if no peaks.
        """
        if not self.peaks:
            return None
        return max(self.peaks, key=lambda p: p.intensity)
    
    def get_molecular_ion_peak(self) -> Optional[Peak]:
        """Get the molecular ion peak (M+).
        
        Searches for the peak with highest m/z that has significant intensity.
        For EI spectra, this is typically near the exact mass.
        
        Returns:
            Molecular ion peak, or None if not found.
        """
        if not self.peaks:
            return None
        # Find peaks with intensity > 1% relative
        significant_peaks = [p for p in self.peaks if p.intensity > 1.0]
        if not significant_peaks:
            return None
        return max(significant_peaks, key=lambda p: p.mz)
    
    def normalize(self, max_intensity: float = 100.0) -> "QCxMS2Result":
        """Return a new result with normalized peak intensities.
        
        Args:
            max_intensity: Target maximum intensity (default 100.0)
            
        Returns:
            New QCxMS2Result with normalized peaks.
        """
        if not self.peaks:
            return self
        
        current_max = max(p.intensity for p in self.peaks)
        if current_max == 0:
            return self
        
        factor = max_intensity / current_max
        normalized_peaks = [
            Peak(p.mz, p.intensity * factor, p.annotation)
            for p in self.peaks
        ]
        
        normalized_fragments = [
            QCxMS2FragmentInfo(
                fragment_path=f.fragment_path,
                mz=f.mz,
                intensity=f.intensity * factor,
                sum_formula=f.sum_formula,
                fragment_type=f.fragment_type,
                reaction_energy=f.reaction_energy,
                barrier=f.barrier,
            )
            for f in self.fragments
        ]
        
        return QCxMS2Result(
            peaks=normalized_peaks,
            fragments=normalized_fragments,
            mode=self.mode,
            input_file=self.input_file,
            molecular_formula=self.molecular_formula,
            ip_ev=self.ip_ev,
            calculation_level=self.calculation_level,
            version=self.version,
            source_dir=self.source_dir,
            metadata=self.metadata.copy(),
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary."""
        return {
            "peaks": [{"mz": p.mz, "intensity": p.intensity, "annotation": p.annotation} for p in self.peaks],
            "fragments": [
                {
                    "fragment_path": f.fragment_path,
                    "mz": f.mz,
                    "intensity": f.intensity,
                    "sum_formula": f.sum_formula,
                    "fragment_type": f.fragment_type,
                    "reaction_energy": f.reaction_energy,
                    "barrier": f.barrier,
                }
                for f in self.fragments
            ],
            "mode": self.mode.value if self.mode else None,
            "input_file": self.input_file,
            "molecular_formula": self.molecular_formula,
            "ip_ev": self.ip_ev,
            "calculation_level": self.calculation_level,
            "version": self.version,
            "source_dir": self.source_dir,
            "metadata": self.metadata,
        }


def parse_qcxms2_spectrum(text: str, normalize_to: float = 100.0) -> List[Peak]:
    """Parse QCxMS2 spectrum from text.
    
    Parses spectrum data in either space-separated or CSV format.
    Each line should contain: m/z intensity [optional annotation]
    
    QCxMS2 outputs intensities normalized to 10000 (NIST convention).
    This function normalizes to the specified value (default 100).
    
    Args:
        text: Spectrum text with m/z and intensity values
        normalize_to: Target maximum intensity for normalization
        
    Returns:
        List of Peak objects sorted by m/z
        
    Examples:
        >>> text = "31.0 10000\\n46.0 1500\\n29.0 3000"
        >>> peaks = parse_qcxms2_spectrum(text)
        >>> len(peaks)
        3
        >>> peaks[0].mz  # Sorted by m/z
        29.0
    """
    peaks = []
    max_intensity = 0.0
    
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        
        # Handle both space-separated and CSV formats
        if "," in line:
            parts = [p.strip() for p in line.split(",")]
        else:
            parts = line.split()
        
        if len(parts) < 2:
            continue
        
        try:
            mz = float(parts[0])
            intensity = float(parts[1])
            annotation = parts[2] if len(parts) > 2 else None
            
            if intensity > max_intensity:
                max_intensity = intensity
            
            peaks.append(Peak(mz=mz, intensity=intensity, annotation=annotation))
        except (ValueError, IndexError) as e:
            logger.debug(f"Skipping unparseable line: {line} ({e})")
            continue
    
    # Normalize intensities
    if max_intensity > 0 and normalize_to > 0:
        factor = normalize_to / max_intensity
        peaks = [Peak(p.mz, p.intensity * factor, p.annotation) for p in peaks]
    
    # Sort by m/z
    return sorted(peaks, key=lambda p: p.mz)


def parse_qcxms2_results(
    filepath: Union[str, Path],
    normalize_to: float = 100.0,
) -> List[Peak]:
    """Parse QCxMS2 peaks.dat or peaks.csv file.
    
    This is the main function for reading QCxMS2 spectrum output.
    
    Args:
        filepath: Path to peaks.dat or peaks.csv file
        normalize_to: Target maximum intensity (default 100.0)
        
    Returns:
        List of Peak objects sorted by m/z
        
    Raises:
        FileNotFoundError: If the file does not exist
        ValueError: If the file cannot be parsed
        
    Example:
        >>> peaks = parse_qcxms2_results("/path/to/calculation/peaks.dat")
        >>> print(f"Found {len(peaks)} peaks")
        >>> print(f"Base peak at m/z = {max(peaks, key=lambda p: p.intensity).mz}")
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"QCxMS2 results file not found: {filepath}")
    
    with open(filepath, "r") as f:
        content = f.read()
    
    peaks = parse_qcxms2_spectrum(content, normalize_to=normalize_to)
    
    if not peaks:
        raise ValueError(f"No valid peaks found in {filepath}")
    
    logger.info(f"Parsed {len(peaks)} peaks from {filepath}")
    return peaks


def parse_qcxms2_allpeaks(
    filepath: Union[str, Path],
    normalize_to: float = 100.0,
) -> List[QCxMS2FragmentInfo]:
    """Parse QCxMS2 allpeaks.dat file with fragment assignments.
    
    The allpeaks.dat file contains fragment information including:
    - Fragment directory path (e.g., "p1f1", "p2f2")
    - Fragment mass (m/z)
    - Relative intensity
    
    Args:
        filepath: Path to allpeaks.dat file
        normalize_to: Target maximum intensity (default 100.0)
        
    Returns:
        List of QCxMS2FragmentInfo objects with fragment data
        
    Raises:
        FileNotFoundError: If the file does not exist
        
    Example:
        >>> fragments = parse_qcxms2_allpeaks("/path/to/calculation/allpeaks.dat")
        >>> for f in fragments:
        ...     print(f"{f.fragment_path}: m/z={f.mz:.1f}, I={f.intensity:.1f}%")
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"QCxMS2 allpeaks file not found: {filepath}")
    
    fragments = []
    max_intensity = 0.0
    
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            
            # Skip header line and empty lines
            if not line or line.startswith("fragment|") or line.startswith("#"):
                continue
            
            # Parse line: fragment_path  mass  intensity
            # Format varies but typically: "  p1f1  31.018538  10000.0"
            parts = line.split()
            
            if len(parts) < 3:
                continue
            
            try:
                fragment_path = parts[0]
                mz = float(parts[1])
                intensity = float(parts[2])
                
                if intensity > max_intensity:
                    max_intensity = intensity
                
                fragments.append(QCxMS2FragmentInfo(
                    fragment_path=fragment_path,
                    mz=mz,
                    intensity=intensity,
                ))
            except (ValueError, IndexError) as e:
                logger.debug(f"Skipping unparseable allpeaks line: {line} ({e})")
                continue
    
    # Normalize intensities
    if max_intensity > 0 and normalize_to > 0:
        factor = normalize_to / max_intensity
        for f in fragments:
            f.intensity *= factor
    
    logger.info(f"Parsed {len(fragments)} fragment entries from {filepath}")
    return fragments


def parse_qcxms2_fragments_file(
    filepath: Union[str, Path],
) -> List[QCxMS2FragmentInfo]:
    """Parse QCxMS2 fragments file with detailed reaction data.
    
    The fragments file contains information about each fragmentation pathway:
    - Directory path
    - Fragment type (fragmentpair/isomer)
    - Reaction energy, DE, barrier, IRC
    - Fragment mass, formula, and relative intensity
    
    Args:
        filepath: Path to fragments file
        
    Returns:
        List of QCxMS2FragmentInfo objects with reaction data
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"QCxMS2 fragments file not found: {filepath}")
    
    fragments = []
    current_parent: Optional[Tuple[str, float, float]] = None  # (path, reaction_energy, barrier)
    
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            
            # Skip header and empty lines
            if not line or line.startswith("Dir"):
                continue
            
            parts = line.split()
            
            if len(parts) < 2:
                continue
            
            try:
                # Check if this is a parent line (fragmentpair/isomer)
                if "fragmentpair" in line or "isomer" in line:
                    # Parent line: path type reaction_e de barrier irc
                    path = parts[0]
                    frag_type = parts[1]
                    reaction_e = float(parts[2]) if len(parts) > 2 else None
                    barrier = float(parts[4]) if len(parts) > 4 else None
                    current_parent = (path, reaction_e, barrier)
                else:
                    # Fragment line: path mass formula intensity
                    path = parts[0]
                    mz = float(parts[1])
                    sum_formula = parts[2] if len(parts) > 2 else None
                    intensity = float(parts[3]) if len(parts) > 3 else 0.0
                    
                    frag = QCxMS2FragmentInfo(
                        fragment_path=path,
                        mz=mz,
                        intensity=intensity,
                        sum_formula=sum_formula,
                    )
                    
                    if current_parent:
                        frag.reaction_energy = current_parent[1]
                        frag.barrier = current_parent[2]
                    
                    fragments.append(frag)
                    
            except (ValueError, IndexError) as e:
                logger.debug(f"Skipping unparseable fragments line: {line} ({e})")
                continue
    
    return fragments


def _parse_qcxms2_output(output_text: str) -> Dict[str, Any]:
    """Parse QCxMS2 stdout output for metadata extraction.
    
    Args:
        output_text: Content of qcxms2.out or stdout
        
    Returns:
        Dictionary with parsed metadata
    """
    metadata: Dict[str, Any] = {}
    
    # Version pattern
    version_match = re.search(r"QCxMS2 version (\S+)", output_text)
    if version_match:
        metadata["version"] = version_match.group(1)
    
    # Mode detection
    if "-cid" in output_text.lower() or "CID" in output_text:
        metadata["mode"] = "CID"
    elif "-dea" in output_text.lower() or "DEA" in output_text:
        metadata["mode"] = "DEA"
    else:
        metadata["mode"] = "EI"
    
    # IP (ionization potential)
    ip_match = re.search(r"IP of input molecule is\s+([\d.]+)\s*eV", output_text)
    if ip_match:
        metadata["ip_ev"] = float(ip_match.group(1))
    
    # Input file
    input_match = re.search(r"input file given:\s*(\S+)", output_text)
    if input_match:
        metadata["input_file"] = input_match.group(1)
    
    # Calculation level
    level_match = re.search(r"Level for geometry optimizations.*?:\s*(\S+)", output_text)
    if level_match:
        metadata["calculation_level"] = level_match.group(1)
    
    # Number of cores
    cores_match = re.search(r"Number of cores used:\s*(\d+)", output_text)
    if cores_match:
        metadata["num_cores"] = int(cores_match.group(1))
    
    # Energy parameters
    eimp_match = re.search(r"eimp0 is:\s*([\d.]+)\s*eV", output_text)
    if eimp_match:
        metadata["eimp_ev"] = float(eimp_match.group(1))
    
    return metadata


def parse_qcxms2_output_dir(
    directory: Union[str, Path],
    normalize_to: float = 100.0,
) -> QCxMS2Result:
    """Parse a complete QCxMS2 calculation directory.
    
    This function reads all relevant output files from a QCxMS2 calculation
    directory and returns a comprehensive QCxMS2Result object.
    
    Files parsed:
    - peaks.dat or peaks.csv: Main spectrum data
    - allpeaks.dat: Fragment assignments
    - qcxms2.out or stdout: Calculation metadata
    - fragments: Detailed fragment information
    
    Args:
        directory: Path to QCxMS2 calculation directory
        normalize_to: Target maximum intensity (default 100.0)
        
    Returns:
        QCxMS2Result with all parsed data
        
    Raises:
        FileNotFoundError: If directory doesn't exist or lacks spectrum files
        
    Example:
        >>> result = parse_qcxms2_output_dir("/path/to/qcxms2_calc/")
        >>> print(f"Calculation mode: {result.mode.value}")
        >>> print(f"Number of peaks: {len(result.peaks)}")
        >>> print(f"Base peak m/z: {result.get_base_peak().mz}")
    """
    directory = Path(directory)
    
    if not directory.exists():
        raise FileNotFoundError(f"QCxMS2 directory not found: {directory}")
    
    if not directory.is_dir():
        raise ValueError(f"Path is not a directory: {directory}")
    
    result = QCxMS2Result(source_dir=str(directory))
    
    # Parse peaks.dat or peaks.csv
    peaks_file = directory / "peaks.dat"
    if not peaks_file.exists():
        peaks_file = directory / "peaks.csv"
    
    if peaks_file.exists():
        result.peaks = parse_qcxms2_results(peaks_file, normalize_to=normalize_to)
    else:
        logger.warning(f"No peaks.dat or peaks.csv found in {directory}")
    
    # Parse allpeaks.dat for fragment info
    allpeaks_file = directory / "allpeaks.dat"
    if allpeaks_file.exists():
        result.fragments = parse_qcxms2_allpeaks(allpeaks_file, normalize_to=normalize_to)
    
    # Parse fragments file for detailed reaction data
    fragments_file = directory / "fragments"
    if fragments_file.exists():
        detailed_fragments = parse_qcxms2_fragments_file(fragments_file)
        # Merge with existing fragment data
        for df in detailed_fragments:
            # Find matching fragment by path
            for f in result.fragments:
                if f.fragment_path == df.fragment_path:
                    f.sum_formula = df.sum_formula or f.sum_formula
                    f.reaction_energy = df.reaction_energy
                    f.barrier = df.barrier
                    break
            else:
                # Not found, add new entry
                result.fragments.append(df)
    
    # Parse output file for metadata
    output_file = directory / "qcxms2.out"
    if not output_file.exists():
        # Try to find any .out file
        out_files = list(directory.glob("*.out"))
        if out_files:
            output_file = out_files[0]
    
    if output_file.exists():
        try:
            with open(output_file, "r") as f:
                output_text = f.read()
            
            metadata = _parse_qcxms2_output(output_text)
            result.metadata.update(metadata)
            
            # Set main attributes from metadata
            if "version" in metadata:
                result.version = metadata["version"]
            if "input_file" in metadata:
                result.input_file = metadata["input_file"]
            if "ip_ev" in metadata:
                result.ip_ev = metadata["ip_ev"]
            if "calculation_level" in metadata:
                result.calculation_level = metadata["calculation_level"]
            if "mode" in metadata:
                mode_str = metadata["mode"]
                if mode_str == "CID":
                    result.mode = IonizationMode.CID
                elif mode_str == "DEA":
                    result.mode = IonizationMode.DEA
                else:
                    result.mode = IonizationMode.EI
                    
        except Exception as e:
            logger.warning(f"Could not parse output file {output_file}: {e}")
    
    logger.info(f"Parsed QCxMS2 calculation from {directory}: "
                f"{len(result.peaks)} peaks, {len(result.fragments)} fragments")
    
    return result
