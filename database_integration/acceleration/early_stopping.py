"""
Early Stopping for QCxMS2 via DreaMS Embedding Similarity

Monitors QCxMS2 calculation progress and stops early when the calculated
spectrum becomes sufficiently similar to a reference (experimental) spectrum.

This can save significant computation time by avoiding unnecessary exploration
of deep fragmentation levels when the main features are already captured.

Approach:
1. Convert partial QCxMS2 output (allpeaks.dat) to MGF
2. Compute DreaMS embedding for partial spectrum
3. Compare to reference spectrum embedding
4. If similarity > threshold, signal early termination

Expected speedup: 2-5x for molecules where Level 1-2 fragmentation is sufficient.

Usage:
    # Command-line monitoring
    python -m database_integration.acceleration.early_stopping \\
        --calc-dir /path/to/qcxms2_calculation \\
        --reference experimental.mgf \\
        --threshold 0.90

    # Python API
    from database_integration.acceleration.early_stopping import EarlyStoppingMonitor
    
    monitor = EarlyStoppingMonitor(reference_spectrum="experimental.mgf")
    while not monitor.should_stop("/path/to/calculation"):
        time.sleep(30)  # Check every 30 seconds
"""

import os
import sys
import time
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union
from datetime import datetime
import json


@dataclass
class ConvergenceStatus:
    """Status of spectrum convergence."""
    
    converged: bool = False
    similarity: float = 0.0
    n_peaks: int = 0
    fragmentation_level: int = 0
    elapsed_time: float = 0.0
    
    # History
    similarity_history: List[float] = None
    
    def __post_init__(self):
        if self.similarity_history is None:
            self.similarity_history = []
    
    @property
    def is_improving(self) -> bool:
        """Check if similarity is still improving."""
        if len(self.similarity_history) < 2:
            return True
        # Check if last 3 values show improvement
        recent = self.similarity_history[-3:]
        return recent[-1] > recent[0] - 0.01  # Allow small fluctuation


class EarlyStoppingMonitor:
    """
    Monitor QCxMS2 calculation and determine when to stop early.
    
    Uses DreaMS embeddings to compare the evolving calculated spectrum
    against a reference (experimental) spectrum.
    
    Example:
        monitor = EarlyStoppingMonitor(
            reference_spectrum="benzene_experimental.mgf",
            similarity_threshold=0.90
        )
        
        # Check periodically
        while calculation_running:
            status = monitor.check("/path/to/qcxms2_output")
            if status.converged:
                print(f"Converged at level {status.fragmentation_level}")
                terminate_calculation()
                break
            time.sleep(60)
    """
    
    def __init__(
        self,
        reference_spectrum: Optional[str] = None,
        reference_embedding: Optional['np.ndarray'] = None,
        similarity_threshold: float = 0.90,
        min_peaks: int = 3,
        min_improvement: float = 0.01,
        patience: int = 3
    ):
        """
        Initialize early stopping monitor.
        
        Args:
            reference_spectrum: Path to reference MGF file
            reference_embedding: Pre-computed reference embedding (1024-dim)
            similarity_threshold: Stop when similarity exceeds this
            min_peaks: Minimum peaks before checking similarity
            min_improvement: Minimum improvement to continue
            patience: Number of checks without improvement before stopping
        """
        self.reference_spectrum = reference_spectrum
        self.similarity_threshold = similarity_threshold
        self.min_peaks = min_peaks
        self.min_improvement = min_improvement
        self.patience = patience
        
        self._reference_embedding = reference_embedding
        self._dreams_available = None
        self._start_time = None
        self._status = ConvergenceStatus()
        self._stagnation_count = 0
    
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
    
    def _get_reference_embedding(self) -> 'np.ndarray':
        """Get or compute reference embedding."""
        if self._reference_embedding is not None:
            return self._reference_embedding
        
        if not self.reference_spectrum:
            raise ValueError("No reference spectrum or embedding provided")
        
        if not self.dreams_available:
            raise ImportError(
                "DreaMS is required. Run: source /mnt/beegfs/software/dreams/setup-dreams.sh"
            )
        
        from dreams.api import dreams_embeddings
        import numpy as np
        
        emb = dreams_embeddings(self.reference_spectrum)
        if emb.ndim == 1:
            emb = emb[np.newaxis, :]
        
        self._reference_embedding = emb
        return emb
    
    def check(self, calc_dir: str) -> ConvergenceStatus:
        """
        Check convergence status of a QCxMS2 calculation.
        
        Args:
            calc_dir: Path to QCxMS2 calculation directory
            
        Returns:
            ConvergenceStatus with current state
        """
        calc_path = Path(calc_dir)
        
        if self._start_time is None:
            self._start_time = time.time()
        
        # Check for allpeaks.dat (contains current spectrum)
        allpeaks_file = calc_path / "allpeaks.dat"
        if not allpeaks_file.exists():
            # Calculation hasn't produced peaks yet
            return self._status
        
        # Parse current peaks
        peaks = self._parse_allpeaks(allpeaks_file)
        
        if len(peaks) < self.min_peaks:
            self._status.n_peaks = len(peaks)
            return self._status
        
        # Determine fragmentation level from directory structure
        level = self._estimate_fragmentation_level(calc_path)
        
        # Convert to MGF and compute similarity
        similarity = self._compute_similarity(peaks)
        
        # Update status
        self._status.n_peaks = len(peaks)
        self._status.similarity = similarity
        self._status.fragmentation_level = level
        self._status.elapsed_time = time.time() - self._start_time
        self._status.similarity_history.append(similarity)
        
        # Check convergence conditions
        if similarity >= self.similarity_threshold:
            self._status.converged = True
            return self._status
        
        # Check for stagnation
        if not self._status.is_improving:
            self._stagnation_count += 1
            if self._stagnation_count >= self.patience:
                self._status.converged = True  # Stop due to stagnation
        else:
            self._stagnation_count = 0
        
        return self._status
    
    def should_stop(self, calc_dir: str) -> bool:
        """
        Simple check if calculation should stop.
        
        Args:
            calc_dir: Path to QCxMS2 calculation directory
            
        Returns:
            True if calculation should stop
        """
        status = self.check(calc_dir)
        return status.converged
    
    def _parse_allpeaks(self, filepath: Path) -> List[Tuple[float, float]]:
        """
        Parse allpeaks.dat file.
        
        Format:
         fragment|fragment mass|relative intensity
         input structure   78.046950     10000.0
         p2   78.046950       181.7
        """
        peaks = []
        peak_dict = {}  # Aggregate by m/z
        
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('fragment|') or '|' in line:
                    continue
                
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        # Format: fragment_name  mass  intensity
                        mz = float(parts[-2])
                        intensity = float(parts[-1])
                        
                        if intensity > 0:
                            # Aggregate peaks at same m/z
                            mz_rounded = round(mz, 2)
                            if mz_rounded in peak_dict:
                                peak_dict[mz_rounded] += intensity
                            else:
                                peak_dict[mz_rounded] = intensity
                    except (ValueError, IndexError):
                        continue
        
        # Convert to list
        peaks = [(mz, intensity) for mz, intensity in peak_dict.items()]
        return peaks
    
    def _estimate_fragmentation_level(self, calc_path: Path) -> int:
        """Estimate current fragmentation level from directory structure."""
        max_level = 0
        
        for item in calc_path.iterdir():
            if item.is_dir() and item.name.startswith('p'):
                # Count depth: p2p3f1p4 would be level 3
                name = item.name
                level = name.count('p') + name.count('f') - 1
                max_level = max(max_level, level)
        
        return max_level
    
    def _compute_similarity(self, peaks: List[Tuple[float, float]]) -> float:
        """
        Compute DreaMS similarity between current peaks and reference.
        """
        if not self.dreams_available:
            # Fallback: simple cosine similarity on peak vectors
            return self._simple_similarity(peaks)
        
        import numpy as np
        from dreams.api import dreams_embeddings
        
        # Convert peaks to MGF
        mgf_content = self._peaks_to_mgf(peaks)
        
        # Write temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mgf', delete=False) as f:
            f.write(mgf_content)
            temp_path = f.name
        
        try:
            # Compute embedding
            calc_emb = dreams_embeddings(temp_path)
            if calc_emb.ndim == 1:
                calc_emb = calc_emb[np.newaxis, :]
            
            # Get reference embedding
            ref_emb = self._get_reference_embedding()
            
            # Compute cosine similarity
            calc_norm = calc_emb / np.linalg.norm(calc_emb)
            ref_norm = ref_emb / np.linalg.norm(ref_emb)
            similarity = float(np.dot(calc_norm.flatten(), ref_norm.flatten()))
            
            return max(0.0, similarity)  # Clamp to [0, 1]
            
        finally:
            os.unlink(temp_path)
    
    def _simple_similarity(self, peaks: List[Tuple[float, float]]) -> float:
        """
        Fallback: Simple peak-based similarity without DreaMS.
        
        Uses weighted Jaccard-like similarity on binned peaks.
        """
        if not self.reference_spectrum:
            return 0.0
        
        # Load reference peaks
        ref_peaks = self._load_mgf_peaks(self.reference_spectrum)
        
        if not ref_peaks:
            return 0.0
        
        # Bin peaks
        bin_width = 1.0
        calc_bins = {}
        ref_bins = {}
        
        for mz, intensity in peaks:
            bin_idx = int(mz / bin_width)
            calc_bins[bin_idx] = calc_bins.get(bin_idx, 0) + intensity
        
        for mz, intensity in ref_peaks:
            bin_idx = int(mz / bin_width)
            ref_bins[bin_idx] = ref_bins.get(bin_idx, 0) + intensity
        
        # Normalize
        calc_max = max(calc_bins.values()) if calc_bins else 1
        ref_max = max(ref_bins.values()) if ref_bins else 1
        
        calc_bins = {k: v/calc_max for k, v in calc_bins.items()}
        ref_bins = {k: v/ref_max for k, v in ref_bins.items()}
        
        # Compute overlap
        all_bins = set(calc_bins.keys()) | set(ref_bins.keys())
        
        numerator = sum(
            min(calc_bins.get(b, 0), ref_bins.get(b, 0))
            for b in all_bins
        )
        denominator = sum(
            max(calc_bins.get(b, 0), ref_bins.get(b, 0))
            for b in all_bins
        )
        
        return numerator / denominator if denominator > 0 else 0.0
    
    def _peaks_to_mgf(self, peaks: List[Tuple[float, float]]) -> str:
        """Convert peaks to MGF format."""
        # Normalize intensities
        max_int = max(p[1] for p in peaks) if peaks else 1
        
        lines = [
            "BEGIN IONS",
            "TITLE=QCxMS2 partial spectrum",
            f"PEPMASS={max(p[0] for p in peaks):.4f}" if peaks else "PEPMASS=0",
            "CHARGE=1+",
            ""
        ]
        
        for mz, intensity in sorted(peaks, key=lambda p: p[0]):
            norm_int = 100.0 * intensity / max_int
            if norm_int >= 0.01:
                lines.append(f"{mz:.4f} {norm_int:.2f}")
        
        lines.append("END IONS")
        return "\n".join(lines)
    
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
                elif in_spectrum and line and not '=' in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            peaks.append((mz, intensity))
                        except ValueError:
                            pass
        
        return peaks
    
    def get_report(self) -> str:
        """Generate a text report of convergence status."""
        lines = [
            "=" * 50,
            "Early Stopping Monitor Report",
            "=" * 50,
            f"Converged: {self._status.converged}",
            f"Similarity: {self._status.similarity:.4f}",
            f"Threshold: {self.similarity_threshold:.4f}",
            f"Peaks: {self._status.n_peaks}",
            f"Fragmentation Level: {self._status.fragmentation_level}",
            f"Elapsed Time: {self._status.elapsed_time:.1f}s",
            "",
            "Similarity History:",
        ]
        
        for i, sim in enumerate(self._status.similarity_history):
            lines.append(f"  Check {i+1}: {sim:.4f}")
        
        lines.append("=" * 50)
        return "\n".join(lines)


def monitor_calculation(
    calc_dir: str,
    reference_spectrum: str,
    threshold: float = 0.90,
    check_interval: int = 60,
    max_time: int = 3600,
    verbose: bool = True
) -> ConvergenceStatus:
    """
    Monitor a running QCxMS2 calculation for early stopping.
    
    Args:
        calc_dir: Path to QCxMS2 calculation directory
        reference_spectrum: Path to reference MGF file
        threshold: Similarity threshold for convergence
        check_interval: Seconds between checks
        max_time: Maximum monitoring time in seconds
        verbose: Print status updates
        
    Returns:
        Final ConvergenceStatus
    """
    monitor = EarlyStoppingMonitor(
        reference_spectrum=reference_spectrum,
        similarity_threshold=threshold
    )
    
    start_time = time.time()
    
    while True:
        status = monitor.check(calc_dir)
        
        if verbose:
            print(f"[{datetime.now().strftime('%H:%M:%S')}] "
                  f"Similarity: {status.similarity:.4f}, "
                  f"Peaks: {status.n_peaks}, "
                  f"Level: {status.fragmentation_level}")
        
        if status.converged:
            if verbose:
                print(f"\nConverged! Final similarity: {status.similarity:.4f}")
            break
        
        if time.time() - start_time > max_time:
            if verbose:
                print(f"\nMax time reached ({max_time}s)")
            break
        
        # Check if calculation finished
        finished_marker = Path(calc_dir) / "QCxMS2_finished"
        if finished_marker.exists():
            if verbose:
                print("\nCalculation finished naturally")
            break
        
        time.sleep(check_interval)
    
    if verbose:
        print(monitor.get_report())
    
    return status


# CLI
def main():
    """Command-line interface for early stopping monitor."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Monitor QCxMS2 calculation for early stopping"
    )
    parser.add_argument(
        "--calc-dir", "-d",
        required=True,
        help="Path to QCxMS2 calculation directory"
    )
    parser.add_argument(
        "--reference", "-r",
        required=True,
        help="Path to reference spectrum (MGF format)"
    )
    parser.add_argument(
        "--threshold", "-t",
        type=float,
        default=0.90,
        help="Similarity threshold for convergence (default: 0.90)"
    )
    parser.add_argument(
        "--interval", "-i",
        type=int,
        default=60,
        help="Check interval in seconds (default: 60)"
    )
    parser.add_argument(
        "--max-time", "-m",
        type=int,
        default=3600,
        help="Maximum monitoring time in seconds (default: 3600)"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress output"
    )
    
    args = parser.parse_args()
    
    status = monitor_calculation(
        calc_dir=args.calc_dir,
        reference_spectrum=args.reference,
        threshold=args.threshold,
        check_interval=args.interval,
        max_time=args.max_time,
        verbose=not args.quiet
    )
    
    # Exit with appropriate code
    sys.exit(0 if status.converged else 1)


if __name__ == "__main__":
    main()
