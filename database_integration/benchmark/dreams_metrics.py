"""
DreaMS-Based Validation Metrics for QCxMS2

Uses DreaMS (Deep Representations Empowering the Annotation of Mass Spectra)
learned embeddings for more meaningful calculated-vs-experimental spectrum comparison.

DreaMS embeddings (1024-dimensional) capture structural and chemical information
that goes beyond simple peak matching, enabling:
- Structure-aware similarity (similar molecules = similar embeddings)
- Robust to peak shifts and intensity differences
- Library search across 201M spectra database

This module integrates with the existing benchmark framework but provides
superior similarity metrics based on learned representations.

Reference:
    Bushuiev et al., "Emergence of molecular structures from repository-scale
    self-supervised learning on tandem mass spectra", Nature Biotechnology (2025)
    DOI: 10.1038/s41587-025-02663-3

Usage:
    from database_integration.benchmark.dreams_metrics import (
        DreaMSSimilarity,
        dreams_embedding_similarity,
        dreams_library_search,
    )
    
    # Compare two spectra
    calc = DreaMSSimilarity()
    similarity = calc.compare("calculated.mgf", "experimental.mgf")
    
    # Batch comparison
    results = calc.compare_batch(
        calculated_files=["mol1.mgf", "mol2.mgf"],
        reference_files=["exp1.mgf", "exp2.mgf"]
    )
    
    # Library search
    matches = calc.library_search("calculated.mgf", library_dir="massbank/")
"""

import os
import sys
import json
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union, Any
import warnings


@dataclass
class DreaMSComparisonResult:
    """Result of a DreaMS-based spectrum comparison."""
    
    # Similarity scores
    embedding_similarity: float = 0.0  # Cosine in embedding space
    
    # Optional: traditional metrics for comparison
    cosine_similarity: Optional[float] = None
    matched_peaks: Optional[int] = None
    
    # Embedding details
    embedding_dim: int = 1024
    
    # Metadata
    calculated_file: Optional[str] = None
    reference_file: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "embedding_similarity": self.embedding_similarity,
            "cosine_similarity": self.cosine_similarity,
            "matched_peaks": self.matched_peaks,
            "embedding_dim": self.embedding_dim,
            "calculated_file": self.calculated_file,
            "reference_file": self.reference_file,
        }
    
    def __str__(self) -> str:
        return f"DreaMS: {self.embedding_similarity:.4f}"


@dataclass
class LibrarySearchResult:
    """Result of a library search."""
    
    query_file: str
    matches: List[Dict[str, Any]] = field(default_factory=list)
    n_matches: int = 0
    search_database: str = ""
    
    def top_match(self) -> Optional[Dict[str, Any]]:
        """Get best match."""
        return self.matches[0] if self.matches else None
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "query_file": self.query_file,
            "n_matches": self.n_matches,
            "search_database": self.search_database,
            "matches": self.matches,
        }


class DreaMSSimilarity:
    """
    Calculate spectrum similarity using DreaMS embeddings.
    
    DreaMS provides learned 1024-dimensional embeddings that capture
    structural and chemical properties of molecules. Similarity in this
    space is more meaningful than raw peak matching.
    
    Example:
        calc = DreaMSSimilarity()
        
        # Single comparison
        result = calc.compare("qcxms2_output.mgf", "massbank_ref.mgf")
        print(f"DreaMS similarity: {result.embedding_similarity:.4f}")
        
        # Batch comparison
        results = calc.compare_batch(
            calculated_files=["mol1.mgf", "mol2.mgf"],
            reference_files=["exp1.mgf", "exp2.mgf"]
        )
    
    Requirements:
        - DreaMS installation: source /mnt/beegfs/software/dreams/setup-dreams.sh
        - GPU recommended for large-scale processing
    """
    
    def __init__(
        self,
        cache_embeddings: bool = True,
        include_traditional_metrics: bool = True,
        mz_tolerance: float = 0.5,
    ):
        """
        Initialize DreaMS similarity calculator.
        
        Args:
            cache_embeddings: Cache computed embeddings for reuse
            include_traditional_metrics: Also compute cosine similarity
            mz_tolerance: m/z tolerance for peak matching (traditional metrics)
        """
        self.cache_embeddings = cache_embeddings
        self.include_traditional_metrics = include_traditional_metrics
        self.mz_tolerance = mz_tolerance
        
        self._embedding_cache: Dict[str, 'np.ndarray'] = {}
        self._dreams_available = None
        self._dreams_model = None
    
    @property
    def dreams_available(self) -> bool:
        """Check if DreaMS is available."""
        if self._dreams_available is None:
            try:
                from dreams.api import dreams_embeddings
                self._dreams_available = True
            except ImportError:
                self._dreams_available = False
        return self._dreams_available
    
    def _get_embedding(self, mgf_path: str) -> 'np.ndarray':
        """
        Get embedding for a spectrum, using cache if available.
        
        Args:
            mgf_path: Path to MGF file
            
        Returns:
            1024-dimensional numpy array
        """
        import numpy as np
        
        # Check cache
        abs_path = str(Path(mgf_path).resolve())
        if self.cache_embeddings and abs_path in self._embedding_cache:
            return self._embedding_cache[abs_path]
        
        if not self.dreams_available:
            raise ImportError(
                "DreaMS is required but not installed. "
                "Run: source /mnt/beegfs/software/dreams/setup-dreams.sh"
            )
        
        from dreams.api import dreams_embeddings
        
        # Compute embedding
        emb = dreams_embeddings(mgf_path)
        
        # Ensure correct shape
        if emb.ndim == 1:
            emb = emb[np.newaxis, :]
        emb = emb.flatten()
        
        # Cache
        if self.cache_embeddings:
            self._embedding_cache[abs_path] = emb
        
        return emb
    
    def _cosine_similarity(
        self,
        emb1: 'np.ndarray',
        emb2: 'np.ndarray'
    ) -> float:
        """Compute cosine similarity between embeddings."""
        import numpy as np
        
        norm1 = np.linalg.norm(emb1)
        norm2 = np.linalg.norm(emb2)
        
        if norm1 == 0 or norm2 == 0:
            return 0.0
        
        return float(np.dot(emb1, emb2) / (norm1 * norm2))
    
    def compare(
        self,
        calculated_mgf: str,
        reference_mgf: str,
    ) -> DreaMSComparisonResult:
        """
        Compare calculated and reference spectra using DreaMS embeddings.
        
        Args:
            calculated_mgf: Path to calculated spectrum (QCxMS2 output)
            reference_mgf: Path to reference spectrum (experimental)
            
        Returns:
            DreaMSComparisonResult with similarity scores
        """
        result = DreaMSComparisonResult(
            calculated_file=calculated_mgf,
            reference_file=reference_mgf,
        )
        
        try:
            # Get embeddings
            calc_emb = self._get_embedding(calculated_mgf)
            ref_emb = self._get_embedding(reference_mgf)
            
            # Compute embedding similarity
            result.embedding_similarity = self._cosine_similarity(calc_emb, ref_emb)
            
        except Exception as e:
            warnings.warn(f"DreaMS embedding failed: {e}")
            result.embedding_similarity = 0.0
        
        # Optionally include traditional metrics
        if self.include_traditional_metrics:
            try:
                trad = self._traditional_similarity(calculated_mgf, reference_mgf)
                result.cosine_similarity = trad.get("cosine", None)
                result.matched_peaks = trad.get("matched_peaks", None)
            except Exception:
                pass
        
        return result
    
    def compare_batch(
        self,
        calculated_files: List[str],
        reference_files: List[str],
    ) -> List[DreaMSComparisonResult]:
        """
        Compare multiple pairs of spectra.
        
        Args:
            calculated_files: List of calculated spectrum paths
            reference_files: List of reference spectrum paths
            
        Returns:
            List of comparison results
        """
        if len(calculated_files) != len(reference_files):
            raise ValueError("Lists must have same length")
        
        results = []
        for calc, ref in zip(calculated_files, reference_files):
            result = self.compare(calc, ref)
            results.append(result)
        
        return results
    
    def compare_from_peaks(
        self,
        calculated_peaks: List[Tuple[float, float]],
        reference_mgf: str,
        metadata: Optional[Dict[str, str]] = None,
    ) -> DreaMSComparisonResult:
        """
        Compare peak list to reference spectrum.
        
        Useful for comparing QCxMS2 results without saving to file first.
        
        Args:
            calculated_peaks: List of (m/z, intensity) tuples
            reference_mgf: Path to reference spectrum
            metadata: Optional metadata (title, precursor_mz, etc.)
            
        Returns:
            DreaMSComparisonResult
        """
        # Convert peaks to temporary MGF
        mgf_content = self._peaks_to_mgf(calculated_peaks, metadata)
        
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.mgf', delete=False
        ) as f:
            f.write(mgf_content)
            temp_path = f.name
        
        try:
            result = self.compare(temp_path, reference_mgf)
            result.calculated_file = "<in-memory peaks>"
            return result
        finally:
            os.unlink(temp_path)
    
    def library_search(
        self,
        query_mgf: str,
        library_dir: Optional[str] = None,
        library_files: Optional[List[str]] = None,
        top_k: int = 10,
    ) -> LibrarySearchResult:
        """
        Search for similar spectra in a library.
        
        Args:
            query_mgf: Path to query spectrum (MGF)
            library_dir: Directory containing library MGF files
            library_files: Explicit list of library files
            top_k: Number of top matches to return
            
        Returns:
            LibrarySearchResult with ranked matches
        """
        import numpy as np
        
        # Get query embedding
        query_emb = self._get_embedding(query_mgf)
        
        # Collect library files
        if library_files is None and library_dir is not None:
            lib_path = Path(library_dir)
            library_files = list(lib_path.glob("*.mgf"))
        
        if not library_files:
            return LibrarySearchResult(query_file=query_mgf)
        
        # Compute similarities
        similarities = []
        for lib_file in library_files:
            try:
                lib_emb = self._get_embedding(str(lib_file))
                sim = self._cosine_similarity(query_emb, lib_emb)
                similarities.append({
                    "file": str(lib_file),
                    "similarity": sim,
                    "metadata": self._extract_mgf_metadata(str(lib_file)),
                })
            except Exception as e:
                continue
        
        # Sort by similarity
        similarities.sort(key=lambda x: -x["similarity"])
        
        return LibrarySearchResult(
            query_file=query_mgf,
            matches=similarities[:top_k],
            n_matches=len(similarities),
            search_database=library_dir or "custom",
        )
    
    def _traditional_similarity(
        self,
        calc_mgf: str,
        ref_mgf: str
    ) -> Dict[str, Any]:
        """
        Compute traditional peak-matching similarity.
        """
        calc_peaks = self._load_mgf_peaks(calc_mgf)
        ref_peaks = self._load_mgf_peaks(ref_mgf)
        
        if not calc_peaks or not ref_peaks:
            return {"cosine": 0.0, "matched_peaks": 0}
        
        # Match peaks
        matched = 0
        matched_intensities = []
        
        for c_mz, c_int in calc_peaks:
            for r_mz, r_int in ref_peaks:
                if abs(c_mz - r_mz) <= self.mz_tolerance:
                    matched += 1
                    matched_intensities.append((c_int, r_int))
                    break
        
        # Compute cosine on matched peaks
        if matched_intensities:
            import numpy as np
            calc_vec = np.array([x[0] for x in matched_intensities])
            ref_vec = np.array([x[1] for x in matched_intensities])
            
            norm_c = np.linalg.norm(calc_vec)
            norm_r = np.linalg.norm(ref_vec)
            
            if norm_c > 0 and norm_r > 0:
                cosine = float(np.dot(calc_vec, ref_vec) / (norm_c * norm_r))
            else:
                cosine = 0.0
        else:
            cosine = 0.0
        
        return {
            "cosine": cosine,
            "matched_peaks": matched,
            "total_calc_peaks": len(calc_peaks),
            "total_ref_peaks": len(ref_peaks),
        }
    
    def _load_mgf_peaks(self, filepath: str) -> List[Tuple[float, float]]:
        """Load peaks from MGF file."""
        peaks = []
        in_spectrum = False
        
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if line == "BEGIN IONS":
                    in_spectrum = True
                elif line == "END IONS":
                    in_spectrum = False
                elif in_spectrum and line and '=' not in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            peaks.append((mz, intensity))
                        except ValueError:
                            pass
        
        return peaks
    
    def _peaks_to_mgf(
        self,
        peaks: List[Tuple[float, float]],
        metadata: Optional[Dict[str, str]] = None
    ) -> str:
        """Convert peaks to MGF format."""
        meta = metadata or {}
        
        # Normalize intensities
        max_int = max(p[1] for p in peaks) if peaks else 1
        
        lines = [
            "BEGIN IONS",
            f"TITLE={meta.get('title', 'QCxMS2 calculated spectrum')}",
        ]
        
        if 'precursor_mz' in meta:
            lines.append(f"PEPMASS={meta['precursor_mz']}")
        elif peaks:
            lines.append(f"PEPMASS={max(p[0] for p in peaks):.4f}")
        
        if 'smiles' in meta:
            lines.append(f"SMILES={meta['smiles']}")
        
        lines.append(f"CHARGE={meta.get('charge', '1+')}") 
        lines.append("")
        
        for mz, intensity in sorted(peaks, key=lambda p: p[0]):
            norm_int = 100.0 * intensity / max_int
            if norm_int >= 0.01:
                lines.append(f"{mz:.4f} {norm_int:.2f}")
        
        lines.append("END IONS")
        return "\n".join(lines)
    
    def _extract_mgf_metadata(self, filepath: str) -> Dict[str, str]:
        """Extract metadata from MGF file."""
        metadata = {}
        
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if '=' in line and not line.startswith('#'):
                    key, _, value = line.partition('=')
                    metadata[key.strip()] = value.strip()
                elif line == "END IONS":
                    break
        
        return metadata
    
    def clear_cache(self):
        """Clear embedding cache."""
        self._embedding_cache.clear()
    
    def cache_size(self) -> int:
        """Get number of cached embeddings."""
        return len(self._embedding_cache)


# Convenience functions

def dreams_embedding_similarity(
    calculated_mgf: str,
    reference_mgf: str,
) -> float:
    """
    Quick function to compare two spectra using DreaMS.
    
    Args:
        calculated_mgf: Path to calculated spectrum
        reference_mgf: Path to reference spectrum
        
    Returns:
        DreaMS embedding similarity (0-1)
    """
    calc = DreaMSSimilarity(cache_embeddings=False)
    result = calc.compare(calculated_mgf, reference_mgf)
    return result.embedding_similarity


def dreams_library_search(
    query_mgf: str,
    library_dir: str,
    top_k: int = 10,
) -> List[Dict[str, Any]]:
    """
    Search for similar spectra in a library.
    
    Args:
        query_mgf: Path to query spectrum
        library_dir: Directory with library MGF files
        top_k: Number of matches to return
        
    Returns:
        List of matches with similarity scores
    """
    calc = DreaMSSimilarity(cache_embeddings=True)
    result = calc.library_search(query_mgf, library_dir=library_dir, top_k=top_k)
    return result.matches


def compare_qcxms2_to_massbank(
    peaks_dat: str,
    massbank_mgf: str,
    smiles: Optional[str] = None,
) -> DreaMSComparisonResult:
    """
    Compare QCxMS2 peaks.dat output to MassBank reference.
    
    Convenience function that handles conversion and comparison.
    
    Args:
        peaks_dat: Path to QCxMS2 peaks.dat file
        massbank_mgf: Path to MassBank MGF reference
        smiles: Optional SMILES for metadata
        
    Returns:
        DreaMSComparisonResult
    """
    # Import converter
    try:
        from database_integration.converters import peaks_to_mgf
    except ImportError:
        from .peaks_to_mgf import peaks_to_mgf
    
    # Convert to MGF
    with tempfile.NamedTemporaryFile(
        mode='w', suffix='.mgf', delete=False
    ) as f:
        temp_path = f.name
    
    try:
        peaks_to_mgf(peaks_dat, temp_path, smiles=smiles)
        
        calc = DreaMSSimilarity()
        result = calc.compare(temp_path, massbank_mgf)
        result.calculated_file = peaks_dat
        
        return result
    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


@dataclass
class DreaMSValidationReport:
    """Comprehensive validation report using DreaMS."""
    
    molecule_name: str = ""
    smiles: str = ""
    formula: str = ""
    
    # DreaMS metrics
    dreams_similarity: float = 0.0
    
    # Traditional metrics (for comparison)
    cosine_similarity: float = 0.0
    spectral_entropy: float = 0.0
    matched_peaks: int = 0
    total_calc_peaks: int = 0
    total_ref_peaks: int = 0
    
    # Quality assessment
    quality_grade: str = ""  # A, B, C, D, F
    
    # File paths
    calculated_file: str = ""
    reference_file: str = ""
    
    def __post_init__(self):
        """Assign quality grade based on DreaMS similarity."""
        if self.dreams_similarity >= 0.90:
            self.quality_grade = "A"
        elif self.dreams_similarity >= 0.75:
            self.quality_grade = "B"
        elif self.dreams_similarity >= 0.60:
            self.quality_grade = "C"
        elif self.dreams_similarity >= 0.40:
            self.quality_grade = "D"
        else:
            self.quality_grade = "F"
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "molecule_name": self.molecule_name,
            "smiles": self.smiles,
            "formula": self.formula,
            "dreams_similarity": self.dreams_similarity,
            "cosine_similarity": self.cosine_similarity,
            "spectral_entropy": self.spectral_entropy,
            "matched_peaks": self.matched_peaks,
            "total_calc_peaks": self.total_calc_peaks,
            "total_ref_peaks": self.total_ref_peaks,
            "quality_grade": self.quality_grade,
            "calculated_file": self.calculated_file,
            "reference_file": self.reference_file,
        }
    
    def to_markdown(self) -> str:
        """Generate markdown report."""
        return f"""## Validation Report: {self.molecule_name or 'Unknown'}

**Structure**: {self.smiles or 'N/A'}  
**Formula**: {self.formula or 'N/A'}  
**Quality Grade**: **{self.quality_grade}**

### Similarity Metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| DreaMS Embedding | {self.dreams_similarity:.4f} | {'Excellent' if self.dreams_similarity > 0.9 else 'Good' if self.dreams_similarity > 0.7 else 'Fair' if self.dreams_similarity > 0.5 else 'Poor'} |
| Cosine Similarity | {self.cosine_similarity:.4f} | Traditional peak matching |
| Matched Peaks | {self.matched_peaks}/{self.total_ref_peaks} | Coverage |

### Files
- Calculated: `{self.calculated_file}`
- Reference: `{self.reference_file}`
"""


def generate_validation_report(
    calculated_mgf: str,
    reference_mgf: str,
    molecule_name: str = "",
    smiles: str = "",
    formula: str = "",
) -> DreaMSValidationReport:
    """
    Generate a comprehensive validation report.
    
    Args:
        calculated_mgf: Path to calculated spectrum
        reference_mgf: Path to reference spectrum
        molecule_name: Name of the molecule
        smiles: SMILES string
        formula: Molecular formula
        
    Returns:
        DreaMSValidationReport with all metrics
    """
    calc = DreaMSSimilarity(include_traditional_metrics=True)
    result = calc.compare(calculated_mgf, reference_mgf)
    
    # Get traditional metrics for comparison
    trad = calc._traditional_similarity(calculated_mgf, reference_mgf)
    
    report = DreaMSValidationReport(
        molecule_name=molecule_name,
        smiles=smiles,
        formula=formula,
        dreams_similarity=result.embedding_similarity,
        cosine_similarity=trad.get("cosine", 0.0),
        matched_peaks=trad.get("matched_peaks", 0),
        total_calc_peaks=trad.get("total_calc_peaks", 0),
        total_ref_peaks=trad.get("total_ref_peaks", 0),
        calculated_file=calculated_mgf,
        reference_file=reference_mgf,
    )
    
    return report


# CLI
def main():
    """Command-line interface for DreaMS validation."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Compare spectra using DreaMS embeddings"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Compare command
    compare_parser = subparsers.add_parser(
        "compare",
        help="Compare two spectra"
    )
    compare_parser.add_argument(
        "calculated",
        help="Path to calculated spectrum (MGF)"
    )
    compare_parser.add_argument(
        "reference",
        help="Path to reference spectrum (MGF)"
    )
    compare_parser.add_argument(
        "--name", "-n",
        default="",
        help="Molecule name"
    )
    compare_parser.add_argument(
        "--smiles", "-s",
        default="",
        help="SMILES string"
    )
    compare_parser.add_argument(
        "--output", "-o",
        help="Output file (JSON or MD)"
    )
    
    # Search command
    search_parser = subparsers.add_parser(
        "search",
        help="Search library for similar spectra"
    )
    search_parser.add_argument(
        "query",
        help="Path to query spectrum (MGF)"
    )
    search_parser.add_argument(
        "library",
        help="Path to library directory"
    )
    search_parser.add_argument(
        "--top-k", "-k",
        type=int,
        default=10,
        help="Number of matches to return"
    )
    
    # Batch command
    batch_parser = subparsers.add_parser(
        "batch",
        help="Batch comparison"
    )
    batch_parser.add_argument(
        "--calc-dir",
        required=True,
        help="Directory with calculated spectra"
    )
    batch_parser.add_argument(
        "--ref-dir",
        required=True,
        help="Directory with reference spectra"
    )
    batch_parser.add_argument(
        "--output", "-o",
        default="batch_results.json",
        help="Output file"
    )
    
    args = parser.parse_args()
    
    if args.command == "compare":
        report = generate_validation_report(
            args.calculated,
            args.reference,
            molecule_name=args.name,
            smiles=args.smiles,
        )
        
        if args.output:
            output_path = Path(args.output)
            if output_path.suffix == ".md":
                output_path.write_text(report.to_markdown())
            else:
                output_path.write_text(json.dumps(report.to_dict(), indent=2))
            print(f"Report saved to {args.output}")
        else:
            print(report.to_markdown())
    
    elif args.command == "search":
        matches = dreams_library_search(
            args.query,
            args.library,
            top_k=args.top_k
        )
        
        print(f"\nTop {args.top_k} matches for {args.query}:\n")
        for i, match in enumerate(matches, 1):
            print(f"{i}. {match['file']}")
            print(f"   Similarity: {match['similarity']:.4f}")
            if 'metadata' in match and match['metadata']:
                title = match['metadata'].get('TITLE', 'N/A')
                print(f"   Title: {title}")
            print()
    
    elif args.command == "batch":
        calc_dir = Path(args.calc_dir)
        ref_dir = Path(args.ref_dir)
        
        calc_files = sorted(calc_dir.glob("*.mgf"))
        ref_files = sorted(ref_dir.glob("*.mgf"))
        
        # Match by filename
        results = []
        for calc_file in calc_files:
            ref_file = ref_dir / calc_file.name
            if ref_file.exists():
                result = generate_validation_report(
                    str(calc_file),
                    str(ref_file),
                    molecule_name=calc_file.stem,
                )
                results.append(result.to_dict())
                print(f"{calc_file.name}: DreaMS={result.dreams_similarity:.4f}")
        
        # Save results
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\nResults saved to {args.output}")
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
