#!/usr/bin/env python3
"""
Automated Benchmark Pipeline for QCxMS2

End-to-end pipeline: molecule → QCxMS2 → compare to experimental databases.

This CLI tool automates the validation workflow:
1. Convert SMILES to 3D structure (XYZ)
2. Run QCxMS2 calculation
3. Convert output to MGF
4. Fetch experimental spectra from MassBank/MassSpecGym
5. Compute similarity metrics (traditional + DreaMS)
6. Generate validation report

Usage:
    # Single molecule benchmark
    python run_benchmark.py --smiles "c1ccccc1" --name benzene
    
    # With experimental reference
    python run_benchmark.py --smiles "c1ccccc1" --reference benzene_exp.mgf
    
    # Batch benchmark from CSV
    python run_benchmark.py --batch molecules.csv --output-dir results/
    
    # Generate SLURM job script
    python run_benchmark.py --smiles "c1ccccc1" --generate-slurm benzene_job.slurm

CSV format for batch mode:
    name,smiles,formula,reference
    benzene,c1ccccc1,C6H6,benzene_exp.mgf
    caffeine,Cn1cnc2c1c(=O)n(c(=O)n2C)C,C8H10N4O2,caffeine_exp.mgf
"""

import os
import sys
import json
import csv
import subprocess
import tempfile
import argparse
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple

# Ensure imports work from various locations
try:
    from database_integration.benchmark.dreams_metrics import (
        DreaMSSimilarity,
        generate_validation_report,
        DreaMSValidationReport,
    )
    from database_integration.converters import peaks_to_mgf
except ImportError:
    # Running as standalone script
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    from database_integration.benchmark.dreams_metrics import (
        DreaMSSimilarity,
        generate_validation_report,
        DreaMSValidationReport,
    )
    from database_integration.converters import peaks_to_mgf


@dataclass
class BenchmarkConfig:
    """Configuration for benchmark run."""
    
    # QCxMS2 settings
    mode: str = "ei"  # ei or cid
    ieeatm: float = 0.6
    mlip: bool = True
    nfrag: int = 5  # Maximum fragmentation depth
    
    # SLURM settings
    partition: str = "gpu-all"
    time: str = "04:00:00"
    cpus: int = 16
    memory: str = "32G"
    gpu: bool = True
    
    # Comparison settings
    mz_tolerance: float = 0.5
    use_dreams: bool = True
    
    # Output
    output_dir: str = "./benchmark_results"
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "mode": self.mode,
            "ieeatm": self.ieeatm,
            "mlip": self.mlip,
            "nfrag": self.nfrag,
            "partition": self.partition,
            "time": self.time,
            "cpus": self.cpus,
            "memory": self.memory,
            "gpu": self.gpu,
            "mz_tolerance": self.mz_tolerance,
            "use_dreams": self.use_dreams,
            "output_dir": self.output_dir,
        }


@dataclass
class BenchmarkResult:
    """Result of a benchmark run."""
    
    name: str
    smiles: str
    formula: str = ""
    
    # Status
    success: bool = False
    error_message: str = ""
    
    # Files
    xyz_file: str = ""
    qcxms2_dir: str = ""
    calculated_mgf: str = ""
    reference_mgf: str = ""
    
    # Metrics
    dreams_similarity: float = 0.0
    cosine_similarity: float = 0.0
    matched_peaks: int = 0
    total_peaks: int = 0
    quality_grade: str = ""
    
    # Timing
    qcxms2_time_seconds: float = 0.0
    total_time_seconds: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "smiles": self.smiles,
            "formula": self.formula,
            "success": self.success,
            "error_message": self.error_message,
            "xyz_file": self.xyz_file,
            "qcxms2_dir": self.qcxms2_dir,
            "calculated_mgf": self.calculated_mgf,
            "reference_mgf": self.reference_mgf,
            "dreams_similarity": self.dreams_similarity,
            "cosine_similarity": self.cosine_similarity,
            "matched_peaks": self.matched_peaks,
            "total_peaks": self.total_peaks,
            "quality_grade": self.quality_grade,
            "qcxms2_time_seconds": self.qcxms2_time_seconds,
            "total_time_seconds": self.total_time_seconds,
        }


class BenchmarkPipeline:
    """
    Automated benchmark pipeline for QCxMS2 validation.
    
    Example:
        pipeline = BenchmarkPipeline(config=BenchmarkConfig(mlip=True))
        
        # Single molecule
        result = pipeline.run_single(
            smiles="c1ccccc1",
            name="benzene",
            reference_mgf="benzene_exp.mgf"
        )
        
        # Batch from CSV
        results = pipeline.run_batch("molecules.csv")
        pipeline.generate_report(results)
    """
    
    def __init__(self, config: Optional[BenchmarkConfig] = None):
        """Initialize benchmark pipeline."""
        self.config = config or BenchmarkConfig()
        
        # Ensure output directory exists
        self.output_dir = Path(self.config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # DreaMS calculator
        self._dreams = None
    
    @property
    def dreams(self) -> DreaMSSimilarity:
        """Lazy-load DreaMS calculator."""
        if self._dreams is None:
            self._dreams = DreaMSSimilarity(cache_embeddings=True)
        return self._dreams
    
    def smiles_to_xyz(self, smiles: str, output_file: str) -> bool:
        """
        Convert SMILES to 3D XYZ structure.
        
        Uses Open Babel for structure generation and optimization.
        
        Args:
            smiles: SMILES string
            output_file: Output XYZ file path
            
        Returns:
            True if successful
        """
        try:
            # Try RDKit first (better conformer generation)
            return self._smiles_to_xyz_rdkit(smiles, output_file)
        except ImportError:
            pass
        
        # Fall back to Open Babel
        try:
            cmd = [
                "obabel", "-:", smiles,
                "-oxyz", "-O", output_file,
                "--gen3d", "--best"
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60
            )
            return result.returncode == 0 and Path(output_file).exists()
        except Exception as e:
            print(f"SMILES conversion failed: {e}")
            return False
    
    def _smiles_to_xyz_rdkit(self, smiles: str, output_file: str) -> bool:
        """Convert SMILES to XYZ using RDKit."""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        mol = Chem.AddHs(mol)
        
        # Generate 3D conformer
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Get coordinates
        conf = mol.GetConformer()
        atoms = mol.GetAtoms()
        
        lines = [str(mol.GetNumAtoms()), "Generated by RDKit"]
        
        for atom in atoms:
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            lines.append(f"{symbol}  {pos.x:.6f}  {pos.y:.6f}  {pos.z:.6f}")
        
        with open(output_file, 'w') as f:
            f.write("\n".join(lines))
        
        return True
    
    def run_qcxms2(
        self,
        xyz_file: str,
        work_dir: str,
        timeout: int = 14400
    ) -> Tuple[bool, str]:
        """
        Run QCxMS2 calculation.
        
        Args:
            xyz_file: Input XYZ structure
            work_dir: Working directory
            timeout: Maximum time in seconds
            
        Returns:
            (success, error_message)
        """
        import time
        
        # Build command
        cmd = ["qcxms2", xyz_file]
        
        if self.config.mode == "cid":
            cmd.append("-cid")
        
        if self.config.mlip:
            cmd.append("-mlip")
        
        cmd.extend(["-ieeatm", str(self.config.ieeatm)])
        cmd.extend(["-nfrag", str(self.config.nfrag)])
        
        start_time = time.time()
        
        try:
            result = subprocess.run(
                cmd,
                cwd=work_dir,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            
            elapsed = time.time() - start_time
            
            # Check for peaks.dat
            peaks_file = Path(work_dir) / "peaks.dat"
            if peaks_file.exists():
                return True, ""
            else:
                return False, f"peaks.dat not generated. Exit code: {result.returncode}"
                
        except subprocess.TimeoutExpired:
            return False, f"QCxMS2 timed out after {timeout}s"
        except Exception as e:
            return False, str(e)
    
    def run_single(
        self,
        smiles: str,
        name: str,
        formula: str = "",
        reference_mgf: Optional[str] = None,
        xyz_file: Optional[str] = None,
    ) -> BenchmarkResult:
        """
        Run benchmark for a single molecule.
        
        Args:
            smiles: SMILES string
            name: Molecule name
            formula: Molecular formula (optional)
            reference_mgf: Path to experimental reference (optional)
            xyz_file: Pre-existing XYZ file (optional)
            
        Returns:
            BenchmarkResult
        """
        import time
        
        result = BenchmarkResult(
            name=name,
            smiles=smiles,
            formula=formula,
            reference_mgf=reference_mgf or "",
        )
        
        start_time = time.time()
        
        # Create working directory
        work_dir = self.output_dir / name
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # Step 1: Generate XYZ if not provided
        if xyz_file and Path(xyz_file).exists():
            result.xyz_file = xyz_file
        else:
            xyz_path = work_dir / f"{name}.xyz"
            if not self.smiles_to_xyz(smiles, str(xyz_path)):
                result.error_message = "Failed to convert SMILES to XYZ"
                result.total_time_seconds = time.time() - start_time
                return result
            result.xyz_file = str(xyz_path)
        
        # Step 2: Run QCxMS2
        qcxms2_dir = work_dir / "qcxms2"
        qcxms2_dir.mkdir(exist_ok=True)
        
        # Copy XYZ to QCxMS2 directory
        import shutil
        input_xyz = qcxms2_dir / "input.xyz"
        shutil.copy(result.xyz_file, input_xyz)
        
        result.qcxms2_dir = str(qcxms2_dir)
        
        qcxms2_start = time.time()
        success, error = self.run_qcxms2(str(input_xyz), str(qcxms2_dir))
        result.qcxms2_time_seconds = time.time() - qcxms2_start
        
        if not success:
            result.error_message = f"QCxMS2 failed: {error}"
            result.total_time_seconds = time.time() - start_time
            return result
        
        # Step 3: Convert output to MGF
        peaks_file = qcxms2_dir / "peaks.dat"
        mgf_file = work_dir / f"{name}_calculated.mgf"
        
        try:
            peaks_to_mgf(
                str(peaks_file),
                str(mgf_file),
                title=f"{name} calculated by QCxMS2",
                smiles=smiles,
                formula=formula,
            )
            result.calculated_mgf = str(mgf_file)
        except Exception as e:
            result.error_message = f"MGF conversion failed: {e}"
            result.total_time_seconds = time.time() - start_time
            return result
        
        # Step 4: Compute similarity if reference provided
        if reference_mgf and Path(reference_mgf).exists():
            try:
                if self.config.use_dreams:
                    report = generate_validation_report(
                        str(mgf_file),
                        reference_mgf,
                        molecule_name=name,
                        smiles=smiles,
                        formula=formula,
                    )
                    result.dreams_similarity = report.dreams_similarity
                    result.cosine_similarity = report.cosine_similarity
                    result.matched_peaks = report.matched_peaks
                    result.total_peaks = report.total_ref_peaks
                    result.quality_grade = report.quality_grade
                else:
                    # Traditional metrics only
                    from .metrics import calculate_cosine_similarity
                    # TODO: implement traditional comparison
                    pass
            except Exception as e:
                result.error_message = f"Comparison failed: {e}"
        
        result.success = True
        result.total_time_seconds = time.time() - start_time
        
        return result
    
    def run_batch(
        self,
        csv_file: str,
        limit: Optional[int] = None
    ) -> List[BenchmarkResult]:
        """
        Run batch benchmark from CSV file.
        
        CSV columns: name, smiles, formula (optional), reference (optional)
        
        Args:
            csv_file: Path to CSV file
            limit: Maximum number of molecules to process
            
        Returns:
            List of BenchmarkResult
        """
        results = []
        
        with open(csv_file) as f:
            reader = csv.DictReader(f)
            
            for i, row in enumerate(reader):
                if limit and i >= limit:
                    break
                
                name = row.get('name', f'mol_{i}')
                smiles = row.get('smiles', '')
                formula = row.get('formula', '')
                reference = row.get('reference', '')
                
                if not smiles:
                    continue
                
                print(f"Processing {name}...")
                
                result = self.run_single(
                    smiles=smiles,
                    name=name,
                    formula=formula,
                    reference_mgf=reference if reference else None,
                )
                
                results.append(result)
                
                # Save intermediate results
                self._save_results(results)
        
        return results
    
    def _save_results(self, results: List[BenchmarkResult]):
        """Save results to JSON."""
        output_file = self.output_dir / "benchmark_results.json"
        
        data = {
            "timestamp": datetime.now().isoformat(),
            "config": self.config.to_dict(),
            "results": [r.to_dict() for r in results],
            "summary": self._compute_summary(results),
        }
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    def _compute_summary(self, results: List[BenchmarkResult]) -> Dict[str, Any]:
        """Compute summary statistics."""
        successful = [r for r in results if r.success]
        
        if not successful:
            return {"n_total": len(results), "n_success": 0}
        
        dreams_sims = [r.dreams_similarity for r in successful if r.dreams_similarity > 0]
        cosine_sims = [r.cosine_similarity for r in successful if r.cosine_similarity > 0]
        
        import numpy as np
        
        return {
            "n_total": len(results),
            "n_success": len(successful),
            "n_failed": len(results) - len(successful),
            "dreams_similarity": {
                "mean": float(np.mean(dreams_sims)) if dreams_sims else 0,
                "std": float(np.std(dreams_sims)) if dreams_sims else 0,
                "min": float(np.min(dreams_sims)) if dreams_sims else 0,
                "max": float(np.max(dreams_sims)) if dreams_sims else 0,
            },
            "cosine_similarity": {
                "mean": float(np.mean(cosine_sims)) if cosine_sims else 0,
                "std": float(np.std(cosine_sims)) if cosine_sims else 0,
            },
            "quality_grades": {
                grade: sum(1 for r in successful if r.quality_grade == grade)
                for grade in ["A", "B", "C", "D", "F"]
            },
        }
    
    def generate_report(
        self,
        results: List[BenchmarkResult],
        output_file: Optional[str] = None
    ) -> str:
        """
        Generate markdown benchmark report.
        
        Args:
            results: List of benchmark results
            output_file: Output file path (optional)
            
        Returns:
            Markdown report string
        """
        summary = self._compute_summary(results)
        
        lines = [
            "# QCxMS2 Benchmark Report",
            "",
            f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Total Molecules**: {summary['n_total']}",
            f"**Successful**: {summary['n_success']}",
            f"**Failed**: {summary['n_failed']}",
            "",
            "## Configuration",
            "",
            f"- Mode: {self.config.mode.upper()}",
            f"- IEE/atom: {self.config.ieeatm}",
            f"- MLIP: {'Enabled' if self.config.mlip else 'Disabled'}",
            f"- Max Depth: {self.config.nfrag}",
            "",
            "## Summary Statistics",
            "",
        ]
        
        if summary['n_success'] > 0:
            ds = summary['dreams_similarity']
            lines.extend([
                "### DreaMS Embedding Similarity",
                "",
                f"- Mean: {ds['mean']:.4f} ± {ds['std']:.4f}",
                f"- Range: [{ds['min']:.4f}, {ds['max']:.4f}]",
                "",
                "### Quality Distribution",
                "",
                "| Grade | Count | Percentage |",
                "|-------|-------|------------|",
            ])
            
            grades = summary['quality_grades']
            for grade in ["A", "B", "C", "D", "F"]:
                count = grades.get(grade, 0)
                pct = 100 * count / summary['n_success'] if summary['n_success'] > 0 else 0
                lines.append(f"| {grade} | {count} | {pct:.1f}% |")
        
        lines.extend([
            "",
            "## Individual Results",
            "",
            "| Molecule | SMILES | DreaMS | Cosine | Grade | Time (s) |",
            "|----------|--------|--------|--------|-------|----------|",
        ])
        
        for r in results:
            status = "✓" if r.success else "✗"
            smiles_short = r.smiles[:20] + "..." if len(r.smiles) > 20 else r.smiles
            lines.append(
                f"| {status} {r.name} | `{smiles_short}` | "
                f"{r.dreams_similarity:.3f} | {r.cosine_similarity:.3f} | "
                f"{r.quality_grade or 'N/A'} | {r.qcxms2_time_seconds:.1f} |"
            )
        
        # Failed molecules
        failed = [r for r in results if not r.success]
        if failed:
            lines.extend([
                "",
                "## Failed Molecules",
                "",
            ])
            for r in failed:
                lines.append(f"- **{r.name}**: {r.error_message}")
        
        lines.extend([
            "",
            "---",
            f"*Generated by QCxMS2 Benchmark Pipeline*",
        ])
        
        report = "\n".join(lines)
        
        if output_file:
            Path(output_file).write_text(report)
        
        return report
    
    def generate_slurm_script(
        self,
        smiles: str,
        name: str,
        output_file: str,
        reference_mgf: Optional[str] = None,
    ) -> str:
        """
        Generate SLURM job script for running benchmark.
        
        Args:
            smiles: SMILES string
            name: Molecule name
            output_file: Output script path
            reference_mgf: Optional reference spectrum
            
        Returns:
            SLURM script content
        """
        script = f'''#!/bin/bash
#SBATCH --job-name=qcxms2-{name}
#SBATCH --partition={self.config.partition}
#SBATCH --time={self.config.time}
#SBATCH --cpus-per-task={self.config.cpus}
#SBATCH --mem={self.config.memory}
{"#SBATCH --gres=gpu:1" if self.config.gpu else ""}
#SBATCH --output={name}_%j.out
#SBATCH --error={name}_%j.err

# Load environment
source /mnt/beegfs/software/qcxms2/setup-qcxms2.sh
source /mnt/beegfs/software/xtb-6.7.1/setup-xtb.sh
source /mnt/beegfs/software/orca-6.1.1/setup-orca.sh

# Create working directory
WORKDIR="${{SLURM_SUBMIT_DIR}}/{name}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "Starting QCxMS2 benchmark for {name}"
echo "SMILES: {smiles}"
echo "Date: $(date)"
echo ""

# Generate XYZ from SMILES
obabel -:"{smiles}" -oxyz -O input.xyz --gen3d --best

# Run QCxMS2
qcxms2 input.xyz {"--cid" if self.config.mode == "cid" else ""} \\
    {"-mlip" if self.config.mlip else ""} \\
    -ieeatm {self.config.ieeatm} \\
    -nfrag {self.config.nfrag}

# Convert to MGF
if [ -f "peaks.dat" ]; then
    python3 -c "
import sys
sys.path.insert(0, '/mnt/beegfs/software/qcxms2/src/qcxms2-source')
from database_integration.converters import peaks_to_mgf
peaks_to_mgf('peaks.dat', '{name}_calculated.mgf', smiles='{smiles}')
print('MGF generated successfully')
"
fi

'''
        
        if reference_mgf:
            script += f'''
# Compare to reference
if [ -f "{name}_calculated.mgf" ]; then
    python3 -c "
import sys
sys.path.insert(0, '/mnt/beegfs/software/qcxms2/src/qcxms2-source')
from database_integration.benchmark.dreams_metrics import generate_validation_report
report = generate_validation_report('{name}_calculated.mgf', '{reference_mgf}', molecule_name='{name}', smiles='{smiles}')
print(report.to_markdown())
"
fi
'''
        
        script += f'''
echo ""
echo "Benchmark completed at $(date)"
'''
        
        Path(output_file).write_text(script)
        os.chmod(output_file, 0o755)
        
        return script


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Automated QCxMS2 benchmark pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Input options
    input_group = parser.add_argument_group("Input")
    input_group.add_argument(
        "--smiles", "-s",
        help="SMILES string for single molecule"
    )
    input_group.add_argument(
        "--name", "-n",
        default="molecule",
        help="Molecule name (default: molecule)"
    )
    input_group.add_argument(
        "--formula", "-f",
        default="",
        help="Molecular formula"
    )
    input_group.add_argument(
        "--reference", "-r",
        help="Reference experimental spectrum (MGF)"
    )
    input_group.add_argument(
        "--xyz",
        help="Pre-existing XYZ file"
    )
    input_group.add_argument(
        "--batch", "-b",
        help="CSV file for batch benchmarking"
    )
    input_group.add_argument(
        "--limit",
        type=int,
        help="Limit number of molecules in batch mode"
    )
    
    # QCxMS2 options
    qcxms2_group = parser.add_argument_group("QCxMS2 Settings")
    qcxms2_group.add_argument(
        "--mode", "-m",
        choices=["ei", "cid"],
        default="ei",
        help="Ionization mode (default: ei)"
    )
    qcxms2_group.add_argument(
        "--ieeatm",
        type=float,
        default=0.6,
        help="Internal energy per atom (default: 0.6)"
    )
    qcxms2_group.add_argument(
        "--mlip",
        action="store_true",
        default=True,
        help="Use MLIP acceleration (default: True)"
    )
    qcxms2_group.add_argument(
        "--no-mlip",
        action="store_true",
        help="Disable MLIP acceleration"
    )
    qcxms2_group.add_argument(
        "--nfrag",
        type=int,
        default=5,
        help="Maximum fragmentation depth (default: 5)"
    )
    
    # Output options
    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--output-dir", "-o",
        default="./benchmark_results",
        help="Output directory (default: ./benchmark_results)"
    )
    output_group.add_argument(
        "--generate-slurm",
        metavar="FILE",
        help="Generate SLURM job script instead of running"
    )
    output_group.add_argument(
        "--report",
        metavar="FILE",
        help="Generate markdown report to file"
    )
    
    # SLURM options (for script generation)
    slurm_group = parser.add_argument_group("SLURM Settings")
    slurm_group.add_argument(
        "--partition",
        default="gpu-all",
        help="SLURM partition (default: gpu-all)"
    )
    slurm_group.add_argument(
        "--time",
        default="04:00:00",
        help="SLURM time limit (default: 04:00:00)"
    )
    slurm_group.add_argument(
        "--cpus",
        type=int,
        default=16,
        help="Number of CPUs (default: 16)"
    )
    slurm_group.add_argument(
        "--memory",
        default="32G",
        help="Memory limit (default: 32G)"
    )
    
    args = parser.parse_args()
    
    # Build config
    config = BenchmarkConfig(
        mode=args.mode,
        ieeatm=args.ieeatm,
        mlip=args.mlip and not args.no_mlip,
        nfrag=args.nfrag,
        output_dir=args.output_dir,
        partition=args.partition,
        time=args.time,
        cpus=args.cpus,
        memory=args.memory,
    )
    
    pipeline = BenchmarkPipeline(config=config)
    
    # Generate SLURM script
    if args.generate_slurm:
        if not args.smiles:
            parser.error("--smiles required with --generate-slurm")
        
        pipeline.generate_slurm_script(
            smiles=args.smiles,
            name=args.name,
            output_file=args.generate_slurm,
            reference_mgf=args.reference,
        )
        print(f"SLURM script generated: {args.generate_slurm}")
        return
    
    # Batch mode
    if args.batch:
        results = pipeline.run_batch(args.batch, limit=args.limit)
        
        # Generate report
        report = pipeline.generate_report(
            results,
            output_file=args.report
        )
        
        if not args.report:
            print(report)
        else:
            print(f"Report saved to {args.report}")
        
        return
    
    # Single molecule mode
    if args.smiles:
        result = pipeline.run_single(
            smiles=args.smiles,
            name=args.name,
            formula=args.formula,
            reference_mgf=args.reference,
            xyz_file=args.xyz,
        )
        
        # Print result
        if result.success:
            print(f"\n✓ Benchmark completed for {result.name}")
            print(f"  QCxMS2 time: {result.qcxms2_time_seconds:.1f}s")
            print(f"  Calculated MGF: {result.calculated_mgf}")
            
            if result.dreams_similarity > 0:
                print(f"\n  DreaMS Similarity: {result.dreams_similarity:.4f}")
                print(f"  Cosine Similarity: {result.cosine_similarity:.4f}")
                print(f"  Quality Grade: {result.quality_grade}")
        else:
            print(f"\n✗ Benchmark failed for {result.name}")
            print(f"  Error: {result.error_message}")
        
        return
    
    # No input provided
    parser.print_help()


if __name__ == "__main__":
    main()
