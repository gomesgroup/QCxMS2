"""
Fragment Prioritization using DreaMS

Uses DreaMS embeddings to predict which fragments are likely to appear in a
mass spectrum, allowing QCxMS2 to prioritize exploration of relevant fragmentation
pathways.

Two main approaches:
1. Spectrum-guided: Query DreaMS with experimental/reference spectrum to find 
   similar spectra, extract common fragment peaks
2. Structure-guided: Use structural fingerprints to find similar molecules,
   then extract their common fragmentation patterns

Expected speedup: 5-10x by focusing on relevant fragments first.
"""

import os
import json
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union
from collections import Counter


@dataclass
class PrioritizedFragment:
    """A fragment predicted to be important based on DreaMS analysis."""
    
    mz: float
    priority_score: float  # 0-1, higher = more likely to appear
    frequency: int = 1  # How many similar spectra contain this fragment
    
    # Optional annotations
    formula: Optional[str] = None
    annotation: Optional[str] = None
    source_spectra: List[int] = field(default_factory=list)
    
    def __repr__(self) -> str:
        formula_str = f", {self.formula}" if self.formula else ""
        return f"PrioritizedFragment(m/z={self.mz:.4f}{formula_str}, score={self.priority_score:.3f})"


class FragmentPrioritizer:
    """
    Use DreaMS to predict expected fragments for a molecule.
    
    This class provides methods to prioritize which fragmentation pathways
    QCxMS2 should explore first, based on:
    - Experimental reference spectra (if available)
    - Structural similarity to known compounds
    - Common fragmentation patterns in similar molecules
    
    Example:
        prioritizer = FragmentPrioritizer()
        
        # From experimental spectrum
        fragments = prioritizer.from_spectrum("experimental.mgf", top_k=20)
        
        # From SMILES structure
        fragments = prioritizer.from_smiles("c1ccccc1", top_k=20)
        
        # Generate CREST constraint file
        prioritizer.write_crest_constraints(fragments, "constraints.inp")
    """
    
    def __init__(
        self,
        dreams_cache_dir: Optional[str] = None,
        min_mz: float = 10.0,
        max_mz: float = 500.0,
        mz_bin_width: float = 0.1,
        min_intensity: float = 0.01
    ):
        """
        Initialize fragment prioritizer.
        
        Args:
            dreams_cache_dir: Directory for DreaMS model cache
            min_mz: Minimum m/z to consider
            max_mz: Maximum m/z to consider
            mz_bin_width: Bin width for m/z grouping
            min_intensity: Minimum relative intensity to consider
        """
        self.dreams_cache_dir = dreams_cache_dir
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.mz_bin_width = mz_bin_width
        self.min_intensity = min_intensity
        
        self._dreams_available = None
        self._massbank_connector = None
    
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
    
    def from_spectrum(
        self,
        spectrum_path: str,
        library_path: Optional[str] = None,
        top_k_similar: int = 50,
        min_frequency: int = 3,
        top_k_fragments: int = 30
    ) -> List[PrioritizedFragment]:
        """
        Prioritize fragments based on similar spectra in a library.
        
        Args:
            spectrum_path: Path to query spectrum (MGF format)
            library_path: Path to reference library (MGF). If None, uses MassBank.
            top_k_similar: Number of similar spectra to analyze
            min_frequency: Minimum occurrences across similar spectra
            top_k_fragments: Maximum fragments to return
            
        Returns:
            List of PrioritizedFragment objects sorted by priority
        """
        if not self.dreams_available:
            raise ImportError(
                "DreaMS is required for spectrum-based prioritization. "
                "Run: source /mnt/beegfs/software/dreams/setup-dreams.sh"
            )
        
        from dreams.api import dreams_embeddings
        import numpy as np
        
        # Load query spectrum
        query_emb = dreams_embeddings(spectrum_path)
        if query_emb.ndim == 1:
            query_emb = query_emb[np.newaxis, :]
        
        # Get library embeddings and spectra
        if library_path is None:
            # Use MassBank connector if no library provided
            library_spectra, library_embs = self._get_massbank_library()
        else:
            library_embs = dreams_embeddings(library_path)
            library_spectra = self._load_mgf_spectra(library_path)
        
        # Compute similarities
        query_norm = query_emb / np.linalg.norm(query_emb, axis=1, keepdims=True)
        library_norm = library_embs / np.linalg.norm(library_embs, axis=1, keepdims=True)
        similarities = (query_norm @ library_norm.T).flatten()
        
        # Get top-k similar spectra indices
        top_indices = np.argsort(similarities)[::-1][:top_k_similar]
        
        # Extract and count fragment peaks
        fragment_counts = self._count_fragments(
            [library_spectra[i] for i in top_indices],
            similarities[top_indices]
        )
        
        # Convert to prioritized fragments
        fragments = self._create_prioritized_fragments(
            fragment_counts,
            min_frequency=min_frequency,
            top_k=top_k_fragments
        )
        
        return fragments
    
    def from_smiles(
        self,
        smiles: str,
        use_massbank: bool = True,
        use_pubchem: bool = False,
        similarity_threshold: float = 0.7,
        top_k_similar: int = 50,
        min_frequency: int = 3,
        top_k_fragments: int = 30
    ) -> List[PrioritizedFragment]:
        """
        Prioritize fragments based on structurally similar molecules.
        
        Uses molecular fingerprints to find similar compounds in databases,
        then extracts common fragmentation patterns from their spectra.
        
        Args:
            smiles: SMILES string of target molecule
            use_massbank: Search MassBank for similar molecules
            use_pubchem: Search PubChem (slower, more comprehensive)
            similarity_threshold: Minimum Tanimoto similarity
            top_k_similar: Number of similar molecules to analyze
            min_frequency: Minimum fragment occurrences
            top_k_fragments: Maximum fragments to return
            
        Returns:
            List of PrioritizedFragment objects sorted by priority
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, DataStructs
        except ImportError:
            raise ImportError(
                "RDKit is required for structure-based prioritization. "
                "Install with: pip install rdkit"
            )
        
        # Generate fingerprint for query molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        query_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        query_mw = Chem.Descriptors.ExactMolWt(mol)
        
        # Collect spectra from similar molecules
        similar_spectra = []
        
        if use_massbank:
            massbank_spectra = self._search_massbank_by_structure(
                smiles, query_fp, query_mw, similarity_threshold, top_k_similar
            )
            similar_spectra.extend(massbank_spectra)
        
        if not similar_spectra:
            # Fallback: use common fragmentation rules
            return self._predict_common_fragments(smiles, mol, top_k_fragments)
        
        # Count fragments across similar spectra
        fragment_counts = self._count_fragments(similar_spectra)
        
        # Convert to prioritized fragments
        fragments = self._create_prioritized_fragments(
            fragment_counts,
            min_frequency=min_frequency,
            top_k=top_k_fragments
        )
        
        return fragments
    
    def from_common_rules(
        self,
        smiles: str,
        top_k: int = 30
    ) -> List[PrioritizedFragment]:
        """
        Predict fragments using common fragmentation rules.
        
        Uses chemical rules to predict likely fragmentation sites:
        - Loss of common neutral molecules (H2O, CO, CO2, NH3, etc.)
        - Cleavage at functional group boundaries
        - Aromatic ring fragmentations
        
        This is a fallback when no similar spectra are available.
        
        Args:
            smiles: SMILES string
            top_k: Maximum fragments to return
            
        Returns:
            List of PrioritizedFragment objects
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors
        except ImportError:
            raise ImportError("RDKit is required. Install with: pip install rdkit")
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        parent_mz = Descriptors.ExactMolWt(mol)
        fragments = []
        
        # Common neutral losses
        neutral_losses = [
            (1.008, "H", 0.9),      # Loss of H
            (2.016, "H2", 0.7),     # Loss of H2
            (15.024, "CH3", 0.8),   # Loss of methyl
            (17.003, "OH", 0.6),    # Loss of OH
            (18.011, "H2O", 0.85),  # Loss of water
            (27.011, "HCN", 0.5),   # Loss of HCN
            (28.006, "CO", 0.75),   # Loss of CO
            (29.003, "CHO", 0.6),   # Loss of CHO
            (30.011, "CH2O", 0.5),  # Loss of formaldehyde
            (43.018, "C2H3O", 0.7), # Loss of acetyl
            (44.010, "CO2", 0.65),  # Loss of CO2
            (45.021, "C2H5O", 0.5), # Loss of ethoxy
        ]
        
        for loss_mass, loss_name, base_priority in neutral_losses:
            product_mz = parent_mz - loss_mass
            if product_mz > self.min_mz:
                fragments.append(PrioritizedFragment(
                    mz=product_mz,
                    priority_score=base_priority,
                    annotation=f"[M-{loss_name}]+"
                ))
        
        # Add parent ion [M]+
        fragments.append(PrioritizedFragment(
            mz=parent_mz,
            priority_score=1.0,
            annotation="[M]+"
        ))
        
        # Detect aromatic rings and add ring fragmentation products
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) == 6:  # Benzene ring
                # Common benzene fragments
                benzene_fragments = [
                    (77.039, "C6H5+", 0.85),   # Phenyl cation
                    (51.023, "C4H3+", 0.7),    # From benzene
                    (39.023, "C3H3+", 0.65),   # Cyclopropenyl
                    (65.039, "C5H5+", 0.6),    # Cyclopentadienyl
                ]
                for mz, formula, priority in benzene_fragments:
                    if mz < parent_mz:
                        fragments.append(PrioritizedFragment(
                            mz=mz,
                            priority_score=priority,
                            formula=formula,
                            annotation="Aromatic ring fragment"
                        ))
                break  # Only process first benzene ring
        
        # Sort by priority and deduplicate
        fragments = self._deduplicate_fragments(fragments)
        fragments.sort(key=lambda f: -f.priority_score)
        
        return fragments[:top_k]
    
    def write_crest_constraints(
        self,
        fragments: List[PrioritizedFragment],
        output_path: str,
        max_fragments: int = 20
    ) -> None:
        """
        Write CREST constraint file for prioritized fragments.
        
        Generates a file that can be used with CREST msreact to focus
        the conformer search on finding specific fragmentation products.
        
        Args:
            fragments: List of prioritized fragments
            output_path: Output file path
            max_fragments: Maximum fragments to include
        """
        selected = fragments[:max_fragments]
        
        with open(output_path, 'w') as f:
            f.write("# CREST msreact constraints generated by QCxMS2 FragmentPrioritizer\n")
            f.write("# Priority-weighted fragment targets\n")
            f.write("$msreact\n")
            
            for frag in selected:
                # Write target mass
                f.write(f"  mstarget={frag.mz:.4f}\n")
                if frag.formula:
                    f.write(f"  # {frag.formula} (priority={frag.priority_score:.3f})\n")
            
            f.write("$end\n")
        
        print(f"Wrote {len(selected)} fragment constraints to {output_path}")
    
    def write_priority_report(
        self,
        fragments: List[PrioritizedFragment],
        output_path: str,
        include_json: bool = True
    ) -> None:
        """
        Write a detailed report of prioritized fragments.
        
        Args:
            fragments: List of prioritized fragments
            output_path: Output file path (markdown)
            include_json: Also write JSON version
        """
        with open(output_path, 'w') as f:
            f.write("# Fragment Priority Report\n\n")
            f.write("Generated by QCxMS2 FragmentPrioritizer using DreaMS\n\n")
            f.write("## Prioritized Fragments\n\n")
            f.write("| Rank | m/z | Formula | Priority | Frequency | Annotation |\n")
            f.write("|------|-----|---------|----------|-----------|------------|\n")
            
            for i, frag in enumerate(fragments, 1):
                formula = frag.formula or "-"
                annotation = frag.annotation or "-"
                f.write(f"| {i} | {frag.mz:.4f} | {formula} | {frag.priority_score:.3f} | {frag.frequency} | {annotation} |\n")
        
        if include_json:
            json_path = output_path.replace('.md', '.json')
            with open(json_path, 'w') as f:
                json.dump([{
                    'mz': frag.mz,
                    'priority_score': frag.priority_score,
                    'frequency': frag.frequency,
                    'formula': frag.formula,
                    'annotation': frag.annotation
                } for frag in fragments], f, indent=2)
    
    def _count_fragments(
        self,
        spectra: List[Dict],
        weights: Optional[List[float]] = None
    ) -> Dict[float, Tuple[float, int, List[int]]]:
        """
        Count fragment occurrences across spectra, grouping by m/z bin.
        
        Returns dict mapping binned m/z to (weighted_count, frequency, source_indices)
        """
        if weights is None:
            weights = [1.0] * len(spectra)
        
        fragment_data = {}  # binned_mz -> (weighted_sum, count, [source_indices])
        
        for spec_idx, (spectrum, weight) in enumerate(zip(spectra, weights)):
            peaks = spectrum.get('peaks', [])
            for mz, intensity in peaks:
                if self.min_mz <= mz <= self.max_mz and intensity >= self.min_intensity:
                    binned_mz = round(mz / self.mz_bin_width) * self.mz_bin_width
                    
                    if binned_mz not in fragment_data:
                        fragment_data[binned_mz] = [0.0, 0, []]
                    
                    fragment_data[binned_mz][0] += weight * intensity
                    fragment_data[binned_mz][1] += 1
                    if spec_idx not in fragment_data[binned_mz][2]:
                        fragment_data[binned_mz][2].append(spec_idx)
        
        return fragment_data
    
    def _create_prioritized_fragments(
        self,
        fragment_counts: Dict[float, Tuple[float, int, List[int]]],
        min_frequency: int = 3,
        top_k: int = 30
    ) -> List[PrioritizedFragment]:
        """Convert fragment counts to prioritized fragment list."""
        fragments = []
        
        # Calculate max weighted count for normalization
        max_weighted = max(d[0] for d in fragment_counts.values()) if fragment_counts else 1.0
        
        for mz, (weighted_sum, count, sources) in fragment_counts.items():
            if count >= min_frequency:
                priority = weighted_sum / max_weighted
                fragments.append(PrioritizedFragment(
                    mz=mz,
                    priority_score=priority,
                    frequency=count,
                    source_spectra=sources
                ))
        
        # Sort by priority
        fragments.sort(key=lambda f: -f.priority_score)
        
        return fragments[:top_k]
    
    def _deduplicate_fragments(
        self,
        fragments: List[PrioritizedFragment],
        tolerance: float = 0.01
    ) -> List[PrioritizedFragment]:
        """Remove duplicate fragments within tolerance."""
        unique = []
        for frag in sorted(fragments, key=lambda f: -f.priority_score):
            is_duplicate = False
            for existing in unique:
                if abs(frag.mz - existing.mz) < tolerance:
                    is_duplicate = True
                    break
            if not is_duplicate:
                unique.append(frag)
        return unique
    
    def _load_mgf_spectra(self, path: str) -> List[Dict]:
        """Load spectra from MGF file."""
        spectra = []
        current_spectrum = {'peaks': []}
        
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line == "BEGIN IONS":
                    current_spectrum = {'peaks': []}
                elif line == "END IONS":
                    spectra.append(current_spectrum)
                elif '=' in line:
                    key, value = line.split('=', 1)
                    current_spectrum[key.lower()] = value
                elif line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            current_spectrum['peaks'].append((mz, intensity))
                        except ValueError:
                            pass
        
        return spectra
    
    def _get_massbank_library(self):
        """Get MassBank library embeddings (cached)."""
        # This would need a pre-computed MassBank embedding cache
        # For now, raise informative error
        raise NotImplementedError(
            "MassBank library search requires pre-computed embeddings. "
            "Please provide a library_path to from_spectrum() or use from_smiles()."
        )
    
    def _search_massbank_by_structure(
        self,
        smiles: str,
        query_fp,
        query_mw: float,
        similarity_threshold: float,
        top_k: int
    ) -> List[Dict]:
        """Search MassBank for structurally similar molecules."""
        # Try to use our MassBank connector
        try:
            from ..connectors.massbank import MassBankConnector
            
            connector = MassBankConnector()
            # Search by exact mass range
            mw_tolerance = 50.0
            results = connector.search_by_mass(
                query_mw - mw_tolerance,
                query_mw + mw_tolerance,
                limit=top_k * 10
            )
            
            # Filter by structural similarity
            from rdkit import Chem
            from rdkit.Chem import AllChem, DataStructs
            
            similar_spectra = []
            for result in results:
                if 'smiles' in result and result['smiles']:
                    mol = Chem.MolFromSmiles(result['smiles'])
                    if mol:
                        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                        similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
                        if similarity >= similarity_threshold:
                            similar_spectra.append({
                                'peaks': result.get('peaks', []),
                                'similarity': similarity,
                                'smiles': result['smiles']
                            })
            
            # Sort by similarity and return top_k
            similar_spectra.sort(key=lambda x: -x['similarity'])
            return similar_spectra[:top_k]
            
        except Exception as e:
            print(f"Warning: MassBank search failed: {e}")
            return []
    
    def _predict_common_fragments(
        self,
        smiles: str,
        mol,
        top_k: int
    ) -> List[PrioritizedFragment]:
        """Fallback: predict fragments using chemical rules."""
        return self.from_common_rules(smiles, top_k)


# Convenience functions

def prioritize_fragments_from_spectrum(
    spectrum_path: str,
    library_path: Optional[str] = None,
    top_k: int = 30
) -> List[PrioritizedFragment]:
    """
    Prioritize fragments based on similar spectra.
    
    Args:
        spectrum_path: Path to query spectrum (MGF)
        library_path: Path to reference library (MGF), optional
        top_k: Number of fragments to return
        
    Returns:
        List of prioritized fragments
    """
    prioritizer = FragmentPrioritizer()
    return prioritizer.from_spectrum(spectrum_path, library_path, top_k_fragments=top_k)


def prioritize_fragments_from_smiles(
    smiles: str,
    top_k: int = 30,
    use_rules: bool = True
) -> List[PrioritizedFragment]:
    """
    Prioritize fragments based on molecular structure.
    
    Args:
        smiles: SMILES string of target molecule
        top_k: Number of fragments to return
        use_rules: Fall back to chemical rules if database search fails
        
    Returns:
        List of prioritized fragments
    """
    prioritizer = FragmentPrioritizer()
    
    try:
        return prioritizer.from_smiles(smiles, top_k_fragments=top_k)
    except Exception as e:
        if use_rules:
            print(f"Database search failed ({e}), using chemical rules")
            return prioritizer.from_common_rules(smiles, top_k)
        raise
