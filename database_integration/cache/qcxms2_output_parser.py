"""
QCxMS2 Output Parser

Parses QCxMS2 calculation output directories to extract fragment data including:
- Fragment structures (XYZ)
- Masses and formulas
- Barriers and energetics
- Pathway information
"""

import os
import re
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict, Tuple


@dataclass
class FragmentData:
    """Data extracted from a single fragment in QCxMS2 output."""
    
    # Identification
    directory: str  # e.g., "p2p7f1"
    parent_directory: str  # e.g., "p2p7" for "p2p7f1"
    
    # Mass and formula
    mz: float
    formula: str  # Normalized to standard order (C, H, N, O, etc.)
    formula_raw: str  # Original QCxMS2 format (e.g., "H5C6")
    
    # Energetics (kcal/mol unless otherwise noted)
    barrier_kcal_mol: Optional[float] = None
    de_kcal_mol: Optional[float] = None
    sumreac_kcal_mol: Optional[float] = None
    
    # Fragment type
    fragment_type: str = "unknown"  # "isomer" or "fragmentpair"
    
    # Structure
    xyz_content: Optional[str] = None
    atom_count: int = 0
    
    # Parent molecule info
    parent_formula: Optional[str] = None
    parent_mz: Optional[float] = None
    parent_smiles: Optional[str] = None
    
    # Calculation metadata
    charge: int = 1
    multiplicity: int = 1
    calculation_level: str = "gfn2"
    
    # Fingerprinting
    structure_hash: Optional[str] = None


@dataclass 
class QCxMS2Output:
    """Complete parsed output from a QCxMS2 calculation."""
    
    # Calculation directory
    base_dir: Path
    
    # Parent molecule
    parent_formula: str = ""
    parent_mz: float = 0.0
    parent_smiles: Optional[str] = None
    
    # All fragments
    fragments: List[FragmentData] = field(default_factory=list)
    
    # Summary statistics
    total_fragments: int = 0
    unique_mz_values: int = 0
    max_fragmentation_level: int = 0
    
    # Calculation settings
    calculation_type: str = "ei"  # "ei" or "cid"
    ieeatm: Optional[float] = None
    mlip_used: bool = False


class QCxMS2OutputParser:
    """
    Parser for QCxMS2 calculation output directories.
    
    Extracts fragment data, barriers, and structures from completed calculations.
    """
    
    def __init__(self, base_dir: str):
        """
        Initialize parser with calculation directory.
        
        Args:
            base_dir: Path to QCxMS2 calculation output directory
        """
        self.base_dir = Path(base_dir)
        if not self.base_dir.exists():
            raise ValueError(f"Directory does not exist: {base_dir}")
        
    def parse(self) -> QCxMS2Output:
        """
        Parse the QCxMS2 output directory.
        
        Returns:
            QCxMS2Output object with all extracted data
        """
        output = QCxMS2Output(base_dir=self.base_dir)
        
        # Parse parent molecule info
        self._parse_parent_info(output)
        
        # Parse allfragments file (contains barriers and energetics)
        allfragments_path = self.base_dir / "allfragments"
        if allfragments_path.exists():
            self._parse_allfragments(output, allfragments_path)
        
        # Parse individual fragment directories for structures
        self._parse_fragment_directories(output)
        
        # Calculate summary statistics
        output.total_fragments = len(output.fragments)
        output.unique_mz_values = len(set(f.mz for f in output.fragments))
        output.max_fragmentation_level = self._calculate_max_level(output)
        
        return output
    
    def _parse_parent_info(self, output: QCxMS2Output) -> None:
        """Parse parent molecule information."""
        # Try to read mass from parent directory
        mass_file = self.base_dir / "mass"
        if mass_file.exists():
            try:
                output.parent_mz = float(mass_file.read_text().strip())
            except ValueError:
                pass
        
        # Try to detect calculation settings
        # Check if MLIP was used (look for orca_trj.xyz in subdirectories)
        for subdir in self.base_dir.iterdir():
            if subdir.is_dir() and (subdir / "orca_trj.xyz").exists():
                # ORCA NEB trajectory exists - likely MLIP-accelerated
                output.mlip_used = True
                break
    
    def _parse_allfragments(self, output: QCxMS2Output, filepath: Path) -> None:
        """
        Parse the allfragments file.
        
        Format (pairs of lines):
        Dir  fragment_type sumreac(kcal/mol)  de(kcal/mol)  barrier(kcal/mol)  irc (cm-1)
         p2 isomer    78.6    78.6    93.4         0.0
         p2     78.112  H6C6         2.1
         
        For fragmentpair type, there are 3 lines:
         p3 fragmentpair    80.2    80.2   152.8         0.0
         p3f1     52.075  H4C4        33.6
         p3f2     26.037  H2C2         0.0
        """
        content = filepath.read_text()
        lines = content.strip().split('\n')
        
        fragments_dict: Dict[str, FragmentData] = {}
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Skip header lines and empty lines
            if line.startswith("Dir") or not line:
                i += 1
                continue
            
            parts = line.split()
            if len(parts) < 2:
                i += 1
                continue
            
            # Check if first part is a directory name (starts with 'p')
            if not parts[0].startswith('p'):
                i += 1
                continue
            
            directory = parts[0]
            
            # Check if this is a pathway header line (has fragment_type)
            if len(parts) >= 6 and parts[1] in ("isomer", "fragmentpair"):
                fragment_type = parts[1]
                try:
                    sumreac = float(parts[2])
                    de = float(parts[3])
                    barrier = float(parts[4])
                except (ValueError, IndexError):
                    i += 1
                    continue
                
                # Read the next line(s) for fragment details
                if fragment_type == "isomer":
                    # Isomer: single fragment with same dir name
                    i += 1
                    if i < len(lines):
                        frag_line = lines[i].strip().split()
                        if len(frag_line) >= 3:
                            frag_dir = frag_line[0]
                            try:
                                mz = float(frag_line[1])
                                formula_raw = frag_line[2]
                                formula = self._normalize_formula(formula_raw)
                                
                                frag = FragmentData(
                                    directory=frag_dir,
                                    parent_directory=self._get_parent_dir(frag_dir),
                                    mz=mz,
                                    formula=formula,
                                    formula_raw=formula_raw,
                                    fragment_type=fragment_type,
                                    barrier_kcal_mol=barrier,
                                    de_kcal_mol=de,
                                    sumreac_kcal_mol=sumreac,
                                    parent_formula=output.parent_formula,
                                    parent_mz=output.parent_mz
                                )
                                fragments_dict[frag_dir] = frag
                            except ValueError:
                                pass
                
                elif fragment_type == "fragmentpair":
                    # Fragment pair: two fragments (f1 and f2)
                    for _ in range(2):
                        i += 1
                        if i < len(lines):
                            frag_line = lines[i].strip().split()
                            if len(frag_line) >= 3:
                                frag_dir = frag_line[0]
                                try:
                                    mz = float(frag_line[1])
                                    formula_raw = frag_line[2]
                                    formula = self._normalize_formula(formula_raw)
                                    
                                    frag = FragmentData(
                                        directory=frag_dir,
                                        parent_directory=self._get_parent_dir(frag_dir),
                                        mz=mz,
                                        formula=formula,
                                        formula_raw=formula_raw,
                                        fragment_type=fragment_type,
                                        barrier_kcal_mol=barrier,
                                        de_kcal_mol=de,
                                        sumreac_kcal_mol=sumreac,
                                        parent_formula=output.parent_formula,
                                        parent_mz=output.parent_mz
                                    )
                                    fragments_dict[frag_dir] = frag
                                except ValueError:
                                    pass
            
            i += 1
        
        output.fragments = list(fragments_dict.values())
    
    def _parse_fragment_directories(self, output: QCxMS2Output) -> None:
        """Parse individual fragment directories for XYZ structures."""
        fragments_by_dir = {f.directory: f for f in output.fragments}
        
        for subdir in self.base_dir.iterdir():
            if not subdir.is_dir():
                continue
            
            # Check if this is a fragment directory
            if subdir.name in fragments_by_dir:
                frag = fragments_by_dir[subdir.name]
                
                # Read fragment.xyz
                xyz_file = subdir / "fragment.xyz"
                if xyz_file.exists():
                    xyz_content = xyz_file.read_text()
                    frag.xyz_content = xyz_content
                    frag.atom_count = self._count_atoms(xyz_content)
                    frag.structure_hash = self._hash_structure(xyz_content)
                
                # Read mass file for more precise value
                mass_file = subdir / "mass"
                if mass_file.exists():
                    try:
                        frag.mz = float(mass_file.read_text().strip())
                    except ValueError:
                        pass
                
                # Read charge
                chrg_file = subdir / ".CHRG"
                if chrg_file.exists():
                    try:
                        frag.charge = int(chrg_file.read_text().strip())
                    except ValueError:
                        pass
                
                # Read multiplicity
                uhf_file = subdir / ".UHF"
                if uhf_file.exists():
                    try:
                        frag.multiplicity = int(uhf_file.read_text().strip()) + 1
                    except ValueError:
                        pass
    
    def _is_pathway_line(self, parts: List[str]) -> bool:
        """Check if parts represent a pathway header line."""
        if len(parts) < 6:
            return False
        # First part should be directory (p2, p2p7, etc.)
        if not re.match(r'^p\d+', parts[0]):
            return False
        # Second part should be fragment_type
        if parts[1] not in ("isomer", "fragmentpair"):
            return False
        return True
    
    def _looks_like_mass(self, s: str) -> bool:
        """Check if string looks like a mass value."""
        try:
            val = float(s)
            return 1.0 < val < 1000.0
        except ValueError:
            return False
    
    def _normalize_formula(self, formula: str) -> str:
        """
        Normalize molecular formula to standard order.
        
        QCxMS2 uses format like "H5C6" but standard is "C6H5"
        """
        # Parse formula into element counts
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)
        
        elements = {}
        for element, count in matches:
            if element:
                count = int(count) if count else 1
                elements[element] = elements.get(element, 0) + count
        
        # Standard order: C, H, then alphabetical
        result = ""
        for elem in ["C", "H"]:
            if elem in elements:
                count = elements.pop(elem)
                result += elem + (str(count) if count > 1 else "")
        
        for elem in sorted(elements.keys()):
            count = elements[elem]
            result += elem + (str(count) if count > 1 else "")
        
        return result
    
    def _get_parent_dir(self, directory: str) -> str:
        """Get parent directory from fragment directory name."""
        # e.g., "p2p7f1" -> "p2p7", "p2p7" -> "p2"
        if directory.endswith(("f1", "f2")):
            return directory[:-2]
        # Remove last pN segment
        match = re.match(r'^(.*)p\d+$', directory)
        if match:
            return match.group(1) or directory
        return directory
    
    def _count_atoms(self, xyz_content: str) -> int:
        """Count atoms in XYZ file content."""
        lines = xyz_content.strip().split('\n')
        if len(lines) >= 1:
            try:
                return int(lines[0].strip())
            except ValueError:
                pass
        return 0
    
    def _hash_structure(self, xyz_content: str) -> str:
        """Create a hash of the atomic structure for matching."""
        # Parse XYZ and create canonical representation
        lines = xyz_content.strip().split('\n')
        if len(lines) < 3:
            return ""
        
        atoms = []
        for line in lines[2:]:  # Skip atom count and comment
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0]
                try:
                    coords = tuple(round(float(x), 2) for x in parts[1:4])
                    atoms.append((element, coords))
                except ValueError:
                    continue
        
        # Sort atoms for canonical ordering
        atoms.sort(key=lambda x: (x[0], x[1]))
        
        # Create hash
        canonical = str(atoms)
        return hashlib.md5(canonical.encode()).hexdigest()[:16]
    
    def _calculate_max_level(self, output: QCxMS2Output) -> int:
        """Calculate maximum fragmentation level from directory names."""
        max_level = 0
        for frag in output.fragments:
            # Count 'p' occurrences minus 1 (first p is level 0)
            level = frag.directory.count('p') - 1
            # Count 'f' occurrences (each f indicates fragmentation)
            level += frag.directory.count('f')
            max_level = max(max_level, level)
        return max_level


def parse_qcxms2_output(directory: str) -> QCxMS2Output:
    """
    Convenience function to parse QCxMS2 output.
    
    Args:
        directory: Path to QCxMS2 calculation directory
        
    Returns:
        QCxMS2Output object with all extracted data
    """
    parser = QCxMS2OutputParser(directory)
    return parser.parse()
