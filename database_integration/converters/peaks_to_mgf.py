"""
QCxMS2 peaks.dat to MGF Format Converter

Converts QCxMS2 output (peaks.dat) to the MGF (Mascot Generic Format) used by
DreaMS, GNPS, and other mass spectrometry tools.

MGF Format Example:
    BEGIN IONS
    TITLE=QCxMS2 calculated spectrum for benzene
    PEPMASS=78.0470
    CHARGE=1+
    IONMODE=Positive
    SOURCE=QCxMS2
    78.0470 100.0
    77.0391 21.3
    52.0313 16.1
    51.0235 16.3
    END IONS

Usage:
    # Single file conversion
    peaks_to_mgf("peaks.dat", "spectrum.mgf", title="benzene")
    
    # Batch conversion
    peaks_to_mgf_batch(
        ["/path/to/calc1/peaks.dat", "/path/to/calc2/peaks.dat"],
        "combined.mgf"
    )
    
    # From Python
    from database_integration.converters import peaks_to_mgf
    peaks_to_mgf("peaks.dat", "output.mgf")
"""

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union
from datetime import datetime


@dataclass
class SpectrumMetadata:
    """Metadata for an MGF spectrum entry."""
    
    title: str = "QCxMS2 calculated spectrum"
    pepmass: Optional[float] = None  # Precursor m/z
    charge: int = 1
    ionmode: str = "Positive"
    source: str = "QCxMS2"
    
    # Optional fields
    name: Optional[str] = None  # Compound name
    formula: Optional[str] = None  # Molecular formula
    smiles: Optional[str] = None  # SMILES string
    inchi: Optional[str] = None  # InChI string
    inchikey: Optional[str] = None  # InChI key
    
    # Calculation metadata
    calculation_type: str = "EI"  # "EI" or "CID"
    ieeatm: Optional[float] = None  # IEE parameter
    calculation_dir: Optional[str] = None  # Source directory
    
    # Custom fields (will be written as-is)
    custom_fields: Dict[str, str] = field(default_factory=dict)


@dataclass
class Peak:
    """A single mass spectrum peak."""
    mz: float
    intensity: float
    annotation: Optional[str] = None


class PeaksToMGFConverter:
    """
    Converts QCxMS2 peaks.dat files to MGF format.
    
    Handles the specific format of QCxMS2 output and generates
    MGF files compatible with DreaMS and other tools.
    """
    
    def __init__(
        self,
        normalize_intensity: bool = True,
        max_intensity: float = 100.0,
        min_intensity: float = 0.01,
        precision_mz: int = 4,
        precision_intensity: int = 2
    ):
        """
        Initialize converter.
        
        Args:
            normalize_intensity: Normalize to base peak = 100
            max_intensity: Maximum intensity value after normalization
            min_intensity: Filter peaks below this intensity
            precision_mz: Decimal places for m/z values
            precision_intensity: Decimal places for intensity values
        """
        self.normalize_intensity = normalize_intensity
        self.max_intensity = max_intensity
        self.min_intensity = min_intensity
        self.precision_mz = precision_mz
        self.precision_intensity = precision_intensity
    
    def parse_peaks_dat(self, filepath: str) -> List[Peak]:
        """
        Parse QCxMS2 peaks.dat file.
        
        Format is two columns: m/z and intensity (can be in scientific notation)
        
        Example:
            78.046949999999981        10000.000000000000
            52.031299999999987        5512.5586768988987
            
        Args:
            filepath: Path to peaks.dat file
            
        Returns:
            List of Peak objects
        """
        peaks = []
        
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        
                        # Skip peaks with negligible intensity
                        if intensity > 0:
                            peaks.append(Peak(mz=mz, intensity=intensity))
                    except ValueError:
                        continue
        
        return peaks
    
    def parse_peaks_csv(self, filepath: str) -> List[Peak]:
        """
        Parse QCxMS2 peaks.csv file (alternative format).
        
        Args:
            filepath: Path to peaks.csv file
            
        Returns:
            List of Peak objects
        """
        peaks = []
        
        with open(filepath) as f:
            header = f.readline()  # Skip header
            for line in f:
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        annotation = parts[2] if len(parts) > 2 else None
                        
                        if intensity > 0:
                            peaks.append(Peak(mz=mz, intensity=intensity, annotation=annotation))
                    except ValueError:
                        continue
        
        return peaks
    
    def normalize_peaks(self, peaks: List[Peak]) -> List[Peak]:
        """
        Normalize peak intensities to base peak = max_intensity.
        
        Args:
            peaks: List of peaks to normalize
            
        Returns:
            New list with normalized intensities
        """
        if not peaks:
            return []
        
        max_int = max(p.intensity for p in peaks)
        if max_int <= 0:
            return peaks
        
        scale_factor = self.max_intensity / max_int
        
        normalized = []
        for peak in peaks:
            new_intensity = peak.intensity * scale_factor
            if new_intensity >= self.min_intensity:
                normalized.append(Peak(
                    mz=peak.mz,
                    intensity=new_intensity,
                    annotation=peak.annotation
                ))
        
        return normalized
    
    def format_mgf(
        self,
        peaks: List[Peak],
        metadata: SpectrumMetadata
    ) -> str:
        """
        Format peaks and metadata as MGF text.
        
        Args:
            peaks: List of peaks
            metadata: Spectrum metadata
            
        Returns:
            MGF formatted string
        """
        lines = ["BEGIN IONS"]
        
        # Add metadata fields
        lines.append(f"TITLE={metadata.title}")
        
        if metadata.pepmass:
            lines.append(f"PEPMASS={metadata.pepmass:.{self.precision_mz}f}")
        elif peaks:
            # Use highest m/z as precursor mass (for EI)
            lines.append(f"PEPMASS={max(p.mz for p in peaks):.{self.precision_mz}f}")
        
        lines.append(f"CHARGE={metadata.charge}+")
        lines.append(f"IONMODE={metadata.ionmode}")
        lines.append(f"SOURCE={metadata.source}")
        
        # Optional standard fields
        if metadata.name:
            lines.append(f"NAME={metadata.name}")
        if metadata.formula:
            lines.append(f"FORMULA={metadata.formula}")
        if metadata.smiles:
            lines.append(f"SMILES={metadata.smiles}")
        if metadata.inchi:
            lines.append(f"INCHI={metadata.inchi}")
        if metadata.inchikey:
            lines.append(f"INCHIKEY={metadata.inchikey}")
        
        # QCxMS2-specific fields
        lines.append(f"CALCTYPE={metadata.calculation_type}")
        if metadata.ieeatm:
            lines.append(f"IEEATM={metadata.ieeatm}")
        if metadata.calculation_dir:
            lines.append(f"CALCDIR={metadata.calculation_dir}")
        
        # Custom fields
        for key, value in metadata.custom_fields.items():
            lines.append(f"{key.upper()}={value}")
        
        lines.append("")  # Empty line before peaks
        
        # Add peaks (sorted by m/z)
        for peak in sorted(peaks, key=lambda p: p.mz):
            mz_str = f"{peak.mz:.{self.precision_mz}f}"
            int_str = f"{peak.intensity:.{self.precision_intensity}f}"
            lines.append(f"{mz_str} {int_str}")
        
        lines.append("END IONS")
        
        return "\n".join(lines)
    
    def convert(
        self,
        input_path: str,
        output_path: str,
        metadata: Optional[SpectrumMetadata] = None
    ) -> int:
        """
        Convert a single peaks.dat file to MGF.
        
        Args:
            input_path: Path to peaks.dat file
            output_path: Path for output MGF file
            metadata: Optional metadata (auto-detected if not provided)
            
        Returns:
            Number of peaks written
        """
        input_path = Path(input_path)
        
        # Parse peaks
        if input_path.suffix == '.csv':
            peaks = self.parse_peaks_csv(str(input_path))
        else:
            peaks = self.parse_peaks_dat(str(input_path))
        
        if not peaks:
            raise ValueError(f"No peaks found in {input_path}")
        
        # Normalize if requested
        if self.normalize_intensity:
            peaks = self.normalize_peaks(peaks)
        
        # Create default metadata if not provided
        if metadata is None:
            metadata = self._auto_metadata(input_path)
        
        # Generate MGF
        mgf_content = self.format_mgf(peaks, metadata)
        
        # Write output
        with open(output_path, 'w') as f:
            f.write(mgf_content)
            f.write("\n")
        
        return len(peaks)
    
    def convert_batch(
        self,
        input_paths: List[str],
        output_path: str,
        metadata_list: Optional[List[SpectrumMetadata]] = None
    ) -> int:
        """
        Convert multiple peaks.dat files to a single MGF file.
        
        Args:
            input_paths: List of paths to peaks.dat files
            output_path: Path for combined output MGF file
            metadata_list: Optional list of metadata (one per input)
            
        Returns:
            Total number of spectra written
        """
        all_mgf = []
        
        for i, input_path in enumerate(input_paths):
            input_path = Path(input_path)
            
            try:
                # Parse peaks
                if input_path.suffix == '.csv':
                    peaks = self.parse_peaks_csv(str(input_path))
                else:
                    peaks = self.parse_peaks_dat(str(input_path))
                
                if not peaks:
                    print(f"Warning: No peaks in {input_path}, skipping")
                    continue
                
                # Normalize
                if self.normalize_intensity:
                    peaks = self.normalize_peaks(peaks)
                
                # Get metadata
                if metadata_list and i < len(metadata_list):
                    metadata = metadata_list[i]
                else:
                    metadata = self._auto_metadata(input_path, index=i)
                
                # Format
                mgf_content = self.format_mgf(peaks, metadata)
                all_mgf.append(mgf_content)
                
            except Exception as e:
                print(f"Warning: Error processing {input_path}: {e}")
                continue
        
        if not all_mgf:
            raise ValueError("No spectra could be converted")
        
        # Write combined output
        with open(output_path, 'w') as f:
            f.write("\n\n".join(all_mgf))
            f.write("\n")
        
        return len(all_mgf)
    
    def _auto_metadata(
        self,
        input_path: Path,
        index: int = 0
    ) -> SpectrumMetadata:
        """
        Auto-generate metadata from file path and directory contents.
        """
        metadata = SpectrumMetadata()
        
        # Try to extract info from directory name
        parent_dir = input_path.parent
        metadata.calculation_dir = str(parent_dir)
        
        # Check for common files that provide metadata
        mass_file = parent_dir / "mass"
        if mass_file.exists():
            try:
                metadata.pepmass = float(mass_file.read_text().strip())
            except:
                pass
        
        # Try to detect calculation type
        if (parent_dir / "cid.inp").exists() or "cid" in str(parent_dir).lower():
            metadata.calculation_type = "CID"
        else:
            metadata.calculation_type = "EI"
        
        # Create title from directory name
        dir_name = parent_dir.name
        metadata.title = f"QCxMS2 {metadata.calculation_type} spectrum - {dir_name}"
        
        # Add timestamp
        metadata.custom_fields["CONVERTED"] = datetime.now().isoformat()
        
        return metadata


# Convenience functions

def peaks_to_mgf(
    input_path: str,
    output_path: str,
    title: Optional[str] = None,
    smiles: Optional[str] = None,
    formula: Optional[str] = None,
    charge: int = 1,
    normalize: bool = True
) -> int:
    """
    Convert QCxMS2 peaks.dat to MGF format.
    
    Args:
        input_path: Path to peaks.dat file
        output_path: Path for output MGF file
        title: Optional spectrum title
        smiles: Optional SMILES string
        formula: Optional molecular formula
        charge: Ion charge (default: 1)
        normalize: Normalize intensities (default: True)
        
    Returns:
        Number of peaks written
        
    Example:
        peaks_to_mgf(
            "benzene_calc/peaks.dat",
            "benzene.mgf",
            title="Benzene EI-MS",
            smiles="c1ccccc1",
            formula="C6H6"
        )
    """
    converter = PeaksToMGFConverter(normalize_intensity=normalize)
    
    metadata = SpectrumMetadata(charge=charge)
    if title:
        metadata.title = title
    if smiles:
        metadata.smiles = smiles
    if formula:
        metadata.formula = formula
    
    return converter.convert(input_path, output_path, metadata)


def peaks_to_mgf_batch(
    input_paths: List[str],
    output_path: str,
    normalize: bool = True
) -> int:
    """
    Convert multiple QCxMS2 peaks.dat files to a single MGF file.
    
    Args:
        input_paths: List of paths to peaks.dat files
        output_path: Path for combined output MGF file
        normalize: Normalize intensities (default: True)
        
    Returns:
        Number of spectra written
        
    Example:
        peaks_to_mgf_batch(
            ["calc1/peaks.dat", "calc2/peaks.dat"],
            "combined.mgf"
        )
    """
    converter = PeaksToMGFConverter(normalize_intensity=normalize)
    return converter.convert_batch(input_paths, output_path)
