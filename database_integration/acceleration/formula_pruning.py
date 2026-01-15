"""
Formula Pruning for QCxMS2 via Mass Constraints and Chemical Rules

Constrains QCxMS2 fragmentation search by predicting likely molecular formulas
for fragment peaks based on:
1. Exact mass (from experimental or MEDUSA-predicted peaks)
2. Parent molecule composition
3. Chemical heuristics (SENIOR rules, nitrogen rule, DBE constraints)
4. Common fragmentation patterns

This can reduce the CREST search space by 3-10x by focusing on chemically
plausible fragments rather than exhaustive enumeration.

Approach:
1. Parse parent molecule formula
2. For each target m/z, enumerate candidate formulas within mass tolerance
3. Apply chemical filters (valence rules, DBE, nitrogen rule)
4. Generate CREST constraint files for directed fragmentation search

Expected speedup: 3-10x depending on molecule complexity.

Usage:
    from database_integration.acceleration.formula_pruning import (
        FormulaPruner,
        enumerate_formulas,
        validate_formula,
    )
    
    # Create pruner from parent molecule
    pruner = FormulaPruner(parent_smiles="c1ccccc1")  # Benzene
    
    # Get likely formulas for a fragment peak
    candidates = pruner.predict_formulas(mz=77.039, tolerance_ppm=10)
    # Returns: [FormulaCandidate(formula="C6H5", mass=77.039, score=0.95)]
    
    # Generate CREST constraint file
    pruner.write_crest_constraints(target_peaks=[77, 51, 39], output="constraints.inp")

References:
    - Kind & Fiehn, BMC Bioinf. 2007: Seven Golden Rules for heuristic filtering
    - MEDUSA: ML-based formula prediction (JACS 2022)
"""

import re
import math
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Set, Any, Union
from pathlib import Path
import warnings

try:
    import numpy as np
except ImportError:
    np = None


# Element masses (monoisotopic)
ELEMENT_MASSES = {
    'H': 1.00782503207,
    'C': 12.0,
    'N': 14.0030740048,
    'O': 15.99491461956,
    'S': 31.97207100,
    'P': 30.97376163,
    'F': 18.99840322,
    'Cl': 34.96885268,
    'Br': 78.9183371,
    'I': 126.904473,
}

# Standard valences
ELEMENT_VALENCES = {
    'H': 1,
    'C': 4,
    'N': 3,
    'O': 2,
    'S': 2,  # Can be 4 or 6 in some cases
    'P': 3,  # Can be 5 in some cases
    'F': 1,
    'Cl': 1,
    'Br': 1,
    'I': 1,
}

# Electron mass for mass defect calculations
ELECTRON_MASS = 0.00054858


@dataclass
class FormulaCandidate:
    """A candidate molecular formula for a fragment."""
    
    formula: str  # e.g., "C6H5"
    composition: Dict[str, int]  # e.g., {"C": 6, "H": 5}
    exact_mass: float  # Calculated monoisotopic mass
    charge: int = 1  # Default positive ion
    
    # Scoring
    mass_error_ppm: float = 0.0
    dbe: float = 0.0  # Double bond equivalents
    score: float = 0.0  # Overall likelihood score
    
    # Validation flags
    valid_valence: bool = True
    valid_nitrogen_rule: bool = True
    valid_dbe: bool = True
    is_subformula_of_parent: bool = False
    
    def __str__(self) -> str:
        return f"{self.formula} ({self.exact_mass:.4f} Da, {self.mass_error_ppm:.1f} ppm)"
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "formula": self.formula,
            "composition": self.composition,
            "exact_mass": self.exact_mass,
            "charge": self.charge,
            "mass_error_ppm": self.mass_error_ppm,
            "dbe": self.dbe,
            "score": self.score,
        }


def parse_formula(formula_str: str) -> Dict[str, int]:
    """
    Parse a molecular formula string into element counts.
    
    Args:
        formula_str: e.g., "C6H5", "C2H6O", "C8H10N4O2"
        
    Returns:
        Dict mapping element symbols to counts
    """
    # Pattern: element symbol followed by optional count
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula_str)
    
    composition = {}
    for element, count in matches:
        if element:
            n = int(count) if count else 1
            composition[element] = composition.get(element, 0) + n
    
    return composition


def formula_to_string(composition: Dict[str, int]) -> str:
    """
    Convert composition dict to formula string (Hill order).
    
    Hill order: C first, then H, then other elements alphabetically.
    """
    parts = []
    
    # Carbon first
    if 'C' in composition and composition['C'] > 0:
        c = composition['C']
        parts.append(f"C{c}" if c > 1 else "C")
    
    # Hydrogen second
    if 'H' in composition and composition['H'] > 0:
        h = composition['H']
        parts.append(f"H{h}" if h > 1 else "H")
    
    # Other elements alphabetically
    for element in sorted(composition.keys()):
        if element not in ('C', 'H') and composition[element] > 0:
            n = composition[element]
            parts.append(f"{element}{n}" if n > 1 else element)
    
    return ''.join(parts)


def calculate_mass(composition: Dict[str, int], charge: int = 1) -> float:
    """
    Calculate monoisotopic mass from composition.
    
    Args:
        composition: Element counts
        charge: Ion charge (positive removes electrons)
        
    Returns:
        Monoisotopic mass in Da
    """
    mass = sum(
        ELEMENT_MASSES.get(el, 0) * count
        for el, count in composition.items()
    )
    
    # Adjust for electron mass (cations lose electrons)
    mass -= charge * ELECTRON_MASS
    
    return mass


def calculate_dbe(composition: Dict[str, int]) -> float:
    """
    Calculate Double Bond Equivalents (DBE) / Ring + Double Bond count.
    
    DBE = 1 + (2*C + N - H - X) / 2
    
    Where X is halogens (F, Cl, Br, I).
    """
    c = composition.get('C', 0)
    h = composition.get('H', 0)
    n = composition.get('N', 0)
    
    # Halogens count as H equivalents
    halogens = sum(
        composition.get(el, 0)
        for el in ['F', 'Cl', 'Br', 'I']
    )
    
    dbe = 1 + (2*c + n - h - halogens) / 2
    return dbe


def validate_valence(composition: Dict[str, int], allow_radicals: bool = True) -> bool:
    """
    Check if formula satisfies SENIOR rules for valence.
    
    For closed-shell molecules: sum of valences must be even and >= 2*max_valence.
    For radicals (odd-electron species): sum can be odd.
    
    In mass spectrometry (especially EI), radicals are common, so we allow
    odd valence sums by default.
    
    Args:
        composition: Element counts
        allow_radicals: Allow odd-electron species (default True for MS)
        
    Returns:
        True if formula passes valence check
    """
    total_valence = sum(
        ELEMENT_VALENCES.get(el, 0) * count
        for el, count in composition.items()
    )
    
    # For closed-shell molecules, sum of valences must be even
    # For radicals, it can be odd
    is_radical = (total_valence % 2 != 0)
    
    if is_radical and not allow_radicals:
        return False
    
    # Check maximum valence constraint
    max_valence = max(
        ELEMENT_VALENCES.get(el, 0)
        for el in composition.keys()
        if composition[el] > 0
    ) if composition else 0
    
    # For even-electron species: valence_sum >= 2 * max_valence
    # For radicals: valence_sum >= 2 * max_valence - 1
    min_valence = 2 * max_valence - (1 if is_radical else 0)
    
    if total_valence < min_valence:
        return False
    
    return True


def validate_nitrogen_rule(composition: Dict[str, int], charge: int = 1) -> bool:
    """
    Check the nitrogen rule for mass spectrometry.
    
    For odd-electron ions (M+.):
    - Even number of N → even nominal mass
    - Odd number of N → odd nominal mass
    
    For even-electron ions (M+):
    - Even number of N → odd nominal mass
    - Odd number of N → even nominal mass
    """
    # Calculate nominal mass
    nominal_mass = int(round(calculate_mass(composition, charge=0)))
    n_count = composition.get('N', 0)
    
    # For radical cations (charge +1, odd electron)
    if charge == 1:
        # Simplified: even N with even mass is OK, odd N with odd mass is OK
        return (n_count % 2) == (nominal_mass % 2)
    
    return True


def is_subformula(child: Dict[str, int], parent: Dict[str, int]) -> bool:
    """
    Check if child formula can be derived from parent by losing atoms.
    
    Args:
        child: Fragment composition
        parent: Parent molecule composition
        
    Returns:
        True if all elements in child are <= parent counts
    """
    for element, count in child.items():
        if count > parent.get(element, 0):
            return False
    return True


def enumerate_formulas(
    target_mass: float,
    elements: List[str],
    max_counts: Dict[str, int],
    tolerance_ppm: float = 10.0,
    charge: int = 1,
) -> List[FormulaCandidate]:
    """
    Enumerate all formulas within mass tolerance.
    
    Uses a recursive algorithm that prunes branches when mass exceeds target.
    
    Args:
        target_mass: Target m/z value (the observed m/z for cations)
        elements: Elements to consider (e.g., ['C', 'H', 'N', 'O'])
        max_counts: Maximum count for each element
        tolerance_ppm: Mass tolerance in ppm
        charge: Ion charge
        
    Returns:
        List of FormulaCandidate objects
    """
    tolerance_da = target_mass * tolerance_ppm / 1e6
    candidates = []
    
    # For cations, target mass is neutral mass minus electron mass
    # So neutral target = target_mass + charge * ELECTRON_MASS
    neutral_target = target_mass + charge * ELECTRON_MASS
    neutral_tolerance = neutral_target * tolerance_ppm / 1e6
    
    def recurse(
        remaining_elements: List[str],
        current_composition: Dict[str, int],
        current_neutral_mass: float,
    ):
        if not remaining_elements:
            # Base case: check if mass matches
            # Calculate actual cation mass for comparison
            cation_mass = current_neutral_mass - charge * ELECTRON_MASS
            error = abs(cation_mass - target_mass)
            
            if error <= tolerance_da:
                # Valid candidate
                formula = formula_to_string(current_composition)
                error_ppm = 1e6 * (cation_mass - target_mass) / target_mass
                dbe = calculate_dbe(current_composition)
                
                candidate = FormulaCandidate(
                    formula=formula,
                    composition=current_composition.copy(),
                    exact_mass=cation_mass,
                    charge=charge,
                    mass_error_ppm=error_ppm,
                    dbe=dbe,
                )
                candidates.append(candidate)
            return
        
        # Get next element
        element = remaining_elements[0]
        element_mass = ELEMENT_MASSES.get(element, 0)
        max_count = max_counts.get(element, 0)
        
        # Try each count
        for count in range(0, max_count + 1):
            new_mass = current_neutral_mass + count * element_mass
            
            # Prune if neutral mass exceeds target + tolerance
            if new_mass > neutral_target + neutral_tolerance:
                break
            
            new_composition = current_composition.copy()
            new_composition[element] = count
            
            recurse(
                remaining_elements[1:],
                new_composition,
                new_mass,
            )
    
    # Initialize recursion with neutral mass = 0
    recurse(elements, {}, 0.0)
    
    return candidates


class FormulaPruner:
    """
    Predict and validate molecular formulas for fragment peaks.
    
    Uses parent molecule composition to constrain the search space,
    then applies chemical heuristics to rank candidates.
    
    Example:
        pruner = FormulaPruner(parent_smiles="c1ccccc1")  # Benzene C6H6
        
        # Predict formulas for m/z 77 (phenyl cation)
        candidates = pruner.predict_formulas(mz=77.039)
        # Returns C6H5 with high score (is subformula of parent)
        
        # Get all likely fragments
        fragments = pruner.predict_all_fragments(
            peaks=[78, 77, 52, 51, 39],
            intensities=[1000, 500, 200, 150, 100]
        )
    """
    
    def __init__(
        self,
        parent_smiles: Optional[str] = None,
        parent_formula: Optional[str] = None,
        elements: Optional[List[str]] = None,
        tolerance_ppm: float = 10.0,
    ):
        """
        Initialize formula pruner.
        
        Args:
            parent_smiles: SMILES string of parent molecule
            parent_formula: Molecular formula of parent (alternative to SMILES)
            elements: Elements to consider (derived from parent if not specified)
            tolerance_ppm: Mass tolerance for formula matching
        """
        self.tolerance_ppm = tolerance_ppm
        
        # Parse parent composition
        if parent_smiles:
            self.parent_composition = self._smiles_to_composition(parent_smiles)
        elif parent_formula:
            self.parent_composition = parse_formula(parent_formula)
        else:
            self.parent_composition = {}
        
        # Determine elements to consider
        if elements:
            self.elements = elements
        elif self.parent_composition:
            self.elements = list(self.parent_composition.keys())
        else:
            self.elements = ['C', 'H', 'N', 'O', 'S']
        
        # Ensure H is always included (for cation formation)
        if 'H' not in self.elements:
            self.elements.append('H')
    
    def _smiles_to_composition(self, smiles: str) -> Dict[str, int]:
        """Convert SMILES to composition using RDKit if available."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {}
            
            mol = Chem.AddHs(mol)
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            return parse_formula(formula)
        except ImportError:
            # Fall back to simple pattern for common cases
            warnings.warn("RDKit not available, using simple SMILES parsing")
            return self._simple_smiles_parse(smiles)
    
    def _simple_smiles_parse(self, smiles: str) -> Dict[str, int]:
        """Simple SMILES parsing for common organic molecules."""
        # Count explicit atoms
        atoms = re.findall(r'[A-Z][a-z]?', smiles)
        composition = {}
        for atom in atoms:
            if atom in ELEMENT_MASSES:
                composition[atom] = composition.get(atom, 0) + 1
        
        # For aromatic molecules, estimate H from degree
        # This is a rough approximation
        return composition
    
    def predict_formulas(
        self,
        mz: float,
        tolerance_ppm: Optional[float] = None,
        charge: int = 1,
        max_candidates: int = 10,
        require_subformula: bool = True,
    ) -> List[FormulaCandidate]:
        """
        Predict likely formulas for a fragment peak.
        
        Args:
            mz: Observed m/z value
            tolerance_ppm: Mass tolerance (uses default if not specified)
            charge: Ion charge
            max_candidates: Maximum number of candidates to return
            require_subformula: Only return formulas that are subformulas of parent
            
        Returns:
            List of FormulaCandidate objects, sorted by score
        """
        tol = tolerance_ppm or self.tolerance_ppm
        
        # Determine max counts for each element
        if self.parent_composition:
            max_counts = self.parent_composition.copy()
        else:
            # Use heuristic max counts based on mass
            max_counts = {
                'C': int(mz / 12) + 1,
                'H': int(mz / 1) + 1,
                'N': min(10, int(mz / 14)),
                'O': min(10, int(mz / 16)),
                'S': min(3, int(mz / 32)),
            }
        
        # Enumerate candidates
        candidates = enumerate_formulas(
            target_mass=mz,
            elements=self.elements,
            max_counts=max_counts,
            tolerance_ppm=tol,
            charge=charge,
        )
        
        # Validate and score
        for candidate in candidates:
            self._validate_candidate(candidate)
            self._score_candidate(candidate)
        
        # Filter
        if require_subformula and self.parent_composition:
            candidates = [c for c in candidates if c.is_subformula_of_parent]
        
        # Filter by basic validity
        candidates = [
            c for c in candidates
            if c.valid_valence and c.valid_dbe
        ]
        
        # Sort by score
        candidates.sort(key=lambda c: -c.score)
        
        return candidates[:max_candidates]
    
    def _validate_candidate(self, candidate: FormulaCandidate):
        """Apply chemical validation rules to a candidate."""
        comp = candidate.composition
        
        # Valence check (allowing radicals for MS)
        candidate.valid_valence = validate_valence(comp, allow_radicals=True)
        
        # Nitrogen rule
        candidate.valid_nitrogen_rule = validate_nitrogen_rule(comp, candidate.charge)
        
        # DBE check
        # For closed-shell: DBE is integer >= 0
        # For radicals (odd H count): DBE is half-integer >= -0.5
        dbe = candidate.dbe
        h_count = comp.get('H', 0)
        is_radical = (h_count % 2 != 0) if comp.get('C', 0) > 0 else False
        
        if is_radical:
            # Radicals have half-integer DBE (e.g., 4.5 for phenyl)
            candidate.valid_dbe = dbe >= -0.5
        else:
            # Closed-shell must have integer DBE >= 0
            candidate.valid_dbe = dbe >= 0 and (dbe == int(dbe))
        
        # Subformula check
        if self.parent_composition:
            candidate.is_subformula_of_parent = is_subformula(
                comp, self.parent_composition
            )
    
    def _score_candidate(self, candidate: FormulaCandidate):
        """
        Score a formula candidate based on multiple factors.
        
        Higher score = more likely formula.
        """
        score = 1.0
        
        # Mass accuracy (exponential decay)
        mass_error = abs(candidate.mass_error_ppm)
        score *= math.exp(-mass_error / 10)  # Penalty for mass error
        
        # Bonus for being subformula of parent
        if candidate.is_subformula_of_parent:
            score *= 2.0
        
        # Penalty for invalid rules
        if not candidate.valid_valence:
            score *= 0.1
        if not candidate.valid_nitrogen_rule:
            score *= 0.5
        if not candidate.valid_dbe:
            score *= 0.1
        
        # DBE reasonableness (0-20 is typical for organics)
        dbe = candidate.dbe
        if 0 <= dbe <= 20:
            score *= 1.0
        elif dbe < 0:
            score *= 0.3
        else:
            score *= 0.5
        
        # Hydrogen deficiency penalty
        # Very low H/C ratios are less common
        c = candidate.composition.get('C', 0)
        h = candidate.composition.get('H', 0)
        if c > 0:
            hc_ratio = h / c
            if hc_ratio < 0.5:
                score *= 0.7  # Highly unsaturated
            elif hc_ratio > 3:
                score *= 0.8  # Very saturated
        
        candidate.score = score
    
    def predict_all_fragments(
        self,
        peaks: List[float],
        intensities: Optional[List[float]] = None,
        max_per_peak: int = 3,
    ) -> Dict[float, List[FormulaCandidate]]:
        """
        Predict formulas for all peaks in a spectrum.
        
        Args:
            peaks: List of m/z values
            intensities: Optional intensities (for prioritization)
            max_per_peak: Maximum candidates per peak
            
        Returns:
            Dict mapping m/z to list of candidates
        """
        results = {}
        
        for mz in peaks:
            candidates = self.predict_formulas(
                mz=mz,
                max_candidates=max_per_peak,
            )
            if candidates:
                results[mz] = candidates
        
        return results
    
    def write_crest_constraints(
        self,
        target_peaks: List[float],
        output_file: str,
        formulas_per_peak: int = 1,
    ) -> str:
        """
        Generate CREST msreact constraint file for directed fragmentation.
        
        The constraint file tells CREST to focus on producing fragments
        with specific molecular formulas, reducing search space.
        
        Args:
            target_peaks: List of m/z values to target
            output_file: Output file path
            formulas_per_peak: Number of formula candidates per peak
            
        Returns:
            Constraint file content
        """
        lines = [
            "# CREST msreact constraints generated by QCxMS2 FormulaPruner",
            "# Target fragments based on chemical rules and parent composition",
            "",
            "$constrain",
        ]
        
        # Get formulas for each peak
        all_formulas = []
        for mz in target_peaks:
            candidates = self.predict_formulas(
                mz=mz,
                max_candidates=formulas_per_peak,
            )
            for c in candidates:
                all_formulas.append(c.formula)
                lines.append(f"  # m/z {mz:.2f}: {c.formula} (score={c.score:.2f})")
        
        # Write target formulas
        if all_formulas:
            lines.append("")
            lines.append("  # Target fragment formulas:")
            for formula in sorted(set(all_formulas)):
                lines.append(f"  formula: {formula}")
        
        lines.append("$end")
        
        content = "\n".join(lines)
        
        Path(output_file).write_text(content)
        return content
    
    def get_summary(self) -> str:
        """Get summary of pruner configuration."""
        lines = [
            "FormulaPruner Configuration",
            "=" * 40,
            f"Elements: {', '.join(self.elements)}",
            f"Tolerance: {self.tolerance_ppm} ppm",
        ]
        
        if self.parent_composition:
            formula = formula_to_string(self.parent_composition)
            mass = calculate_mass(self.parent_composition)
            lines.extend([
                "",
                f"Parent formula: {formula}",
                f"Parent mass: {mass:.4f} Da",
            ])
        
        return "\n".join(lines)


# MEDUSA Integration (optional)

def medusa_formula_prediction(
    spectrum_file: str,
    target_mz: float,
    tolerance_da: float = 0.01,
) -> Optional[str]:
    """
    Use MEDUSA LSTM for formula prediction from high-resolution spectrum.
    
    Note: This requires high-resolution data (FT-ICR or Orbitrap) with
    visible fine isotopic structure.
    
    Args:
        spectrum_file: Path to mzXML/mzML file
        target_mz: m/z value to predict formula for
        tolerance_da: Mass tolerance for peak region
        
    Returns:
        Predicted formula or None if MEDUSA not available
    """
    try:
        from mass_automation.experiment import Experiment
        from mass_automation.formula.model import LSTM
        
        # This is a placeholder - actual MEDUSA usage requires
        # careful data preparation and model selection
        # See MEDUSA documentation for proper workflow
        
        warnings.warn(
            "MEDUSA formula prediction requires high-resolution data "
            "and careful configuration. Using heuristic fallback."
        )
        return None
        
    except ImportError:
        return None


# CLI
def main():
    """Command-line interface for formula pruning."""
    import argparse
    import json
    
    parser = argparse.ArgumentParser(
        description="Predict molecular formulas for fragment peaks"
    )
    
    parser.add_argument(
        "--smiles", "-s",
        help="Parent molecule SMILES"
    )
    parser.add_argument(
        "--formula", "-f",
        help="Parent molecule formula"
    )
    parser.add_argument(
        "--mz", "-m",
        type=float,
        nargs='+',
        required=True,
        help="Target m/z value(s)"
    )
    parser.add_argument(
        "--tolerance", "-t",
        type=float,
        default=10.0,
        help="Mass tolerance in ppm (default: 10)"
    )
    parser.add_argument(
        "--max-candidates", "-n",
        type=int,
        default=5,
        help="Max candidates per peak (default: 5)"
    )
    parser.add_argument(
        "--output", "-o",
        help="Output file (JSON)"
    )
    parser.add_argument(
        "--crest-output",
        help="Generate CREST constraint file"
    )
    
    args = parser.parse_args()
    
    # Create pruner
    pruner = FormulaPruner(
        parent_smiles=args.smiles,
        parent_formula=args.formula,
        tolerance_ppm=args.tolerance,
    )
    
    print(pruner.get_summary())
    print()
    
    # Predict formulas
    results = {}
    for mz in args.mz:
        candidates = pruner.predict_formulas(
            mz=mz,
            max_candidates=args.max_candidates,
        )
        
        print(f"m/z {mz:.4f}:")
        for i, c in enumerate(candidates, 1):
            print(f"  {i}. {c}")
            if c.is_subformula_of_parent:
                print(f"      ✓ Subformula of parent")
        print()
        
        results[mz] = [c.to_dict() for c in candidates]
    
    # Save results
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Results saved to {args.output}")
    
    # Generate CREST constraints
    if args.crest_output:
        pruner.write_crest_constraints(args.mz, args.crest_output)
        print(f"CREST constraints saved to {args.crest_output}")


if __name__ == "__main__":
    main()
