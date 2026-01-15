# QCxMS2 Development Roadmap

> **Branch**: gpg-cluster  
> **Last Updated**: January 15, 2026  
> **Maintainers**: GPG Cluster Team (CMU)

This document tracks planned improvements, integrations, and research directions for QCxMS2 on the GPG cluster.

---

## Table of Contents

1. [Priority Matrix](#priority-matrix)
2. [Phase 1: Quick Wins](#phase-1-quick-wins-this-week)
3. [Phase 2: Validation Workflows](#phase-2-validation-workflows-1-2-weeks)
4. [Phase 3: ML-Accelerated QCxMS2](#phase-3-ml-accelerated-qcxms2-2-4-weeks)
5. [Phase 4: Research Directions](#phase-4-research-directions-future)
6. [Completed Work](#completed-work)
7. [Known Issues](#known-issues)

---

## Priority Matrix

| ID | Task | Impact | Effort | Priority | Status |
|----|------|--------|--------|----------|--------|
| A1 | MLIP for CREST msreact | High | Low | **P0** | In Progress |
| A2 | Fragment m/z cache | Very High | Medium | **P0** | ✅ Complete |
| V1 | QCxMS2-to-MGF converter | Medium | Low | **P1** | Not Started |
| V2 | DreaMS-based validation metric | High | Medium | **P1** | Not Started |
| V3 | Automated benchmark pipeline | High | Medium | **P1** | Partial |
| A3 | Fragment prioritization (DreaMS) | High | Medium | **P1** | Not Started |
| A4 | Formula pruning (MEDUSA LSTM) | Medium | Medium | **P2** | Not Started |
| A5 | Early stopping via embeddings | Medium | Low | **P2** | Not Started |
| A6 | Pre-computed fragment database | Very High | High | **P2** | ✅ Complete (merged with A2) |
| V4 | Multi-molecule benchmark suite | High | High | **P2** | Not Started |
| R1 | ML-predicted NEB barriers | Very High | Very High | **P3** | Research |
| R2 | GNN direct fragmentation | Very High | Very High | **P3** | Research |

**Legend**: P0 = Do Now, P1 = Next Sprint, P2 = This Month, P3 = Future Research

---

## Phase 1: Quick Wins (This Week)

### A1: MLIP Acceleration for CREST msreact
**Status**: In Progress  
**Impact**: 5-10x speedup for conformer/reaction sampling  
**Dependencies**: AIMNet2 server (id-gpu01:8888)

CREST msreact is 30-50% of total QCxMS2 runtime. We already have MLIP for NEB; extend to CREST.

**Implementation**:
```bash
# Already have setup script:
source /mnt/beegfs/software/crest-3.0.2/setup-crest-mlip.sh

# Sets XTB_MLIP_WRAPPER for CREST to use
crest mol.xyz --msreact -xnam $XTB_MLIP_WRAPPER
```

**Tasks**:
- [ ] Benchmark MLIP-CREST vs GFN2-CREST for benzene fragmentation
- [ ] Verify fragment detection rate is equivalent
- [ ] Document in skill file
- [ ] Add to QCxMS2 wrapper script (auto-enable when `-mlip` flag used)

**Files**:
- `/mnt/beegfs/software/crest-3.0.2/setup-crest-mlip.sh`
- `/mnt/beegfs/software/qcxms2/setup-qcxms2.sh` (needs update)

---

### A2: Fragment m/z Cache with medusa-search
**Status**: Not Started  
**Impact**: 50-100x for repeated calculations  
**Dependencies**: MEDUSA installation

Use medusa-search's 4-decimal m/z indexing to cache computed fragments and their NEB barriers.

**Concept**:
```python
# Build index of known fragments
from medusa_search import MedusaIndex

fragment_db = MedusaIndex()
fragment_db.add_batch([
    {"mz": 77.0391, "formula": "C6H5+", "barrier": 3.2, "smiles": "[CH]1=CC=CC=C1"},
    {"mz": 51.0235, "formula": "C4H3+", "barrier": 4.1, "smiles": "..."},
    # ... thousands of pre-computed fragments
])
fragment_db.build_index()
fragment_db.save("fragment_cache.idx")

# In QCxMS2 workflow:
def get_barrier(product_mz, product_structure):
    cached = fragment_db.search(round(product_mz, 4))
    if cached and structure_match(cached.smiles, product_structure):
        return cached.barrier  # Instant!
    else:
        return run_neb(product_structure)  # Only if novel
```

**Tasks**:
- [ ] Design fragment database schema (m/z, formula, barrier, SMILES, source)
- [ ] Create `fragment_cache.py` module
- [ ] Integrate with QCxMS2 NEB step
- [ ] Seed database with fragments from existing calculations
- [ ] Benchmark cache hit rate on test molecules

**Files to create**:
- `database_integration/cache/fragment_cache.py`
- `database_integration/cache/fragment_db.idx`

---

## Phase 2: Validation Workflows (1-2 Weeks)

### V1: QCxMS2-to-MGF Converter
**Status**: Not Started  
**Impact**: Enables DreaMS integration  
**Dependencies**: None

Convert QCxMS2 `peaks.dat` output to MGF format for DreaMS processing.

**Implementation**:
```python
# database_integration/converters/peaks_to_mgf.py

def peaks_to_mgf(peaks_dat: str, output_mgf: str, metadata: dict = None):
    """
    Convert QCxMS2 peaks.dat to MGF format.
    
    Args:
        peaks_dat: Path to QCxMS2 peaks.dat file
        output_mgf: Output MGF file path
        metadata: Optional dict with TITLE, PEPMASS, etc.
    """
    peaks = parse_peaks_dat(peaks_dat)
    
    with open(output_mgf, 'w') as f:
        f.write("BEGIN IONS\n")
        f.write(f"TITLE={metadata.get('title', 'QCxMS2 Calculated')}\n")
        f.write(f"PEPMASS={metadata.get('precursor_mz', peaks[0]['mz'])}\n")
        f.write(f"CHARGE={metadata.get('charge', '1+')}\n")
        for peak in peaks:
            f.write(f"{peak['mz']:.4f} {peak['intensity']:.1f}\n")
        f.write("END IONS\n")
```

**Tasks**:
- [ ] Create `converters/peaks_to_mgf.py`
- [ ] Add batch conversion support
- [ ] Include metadata extraction from QCxMS2 logs
- [ ] Add CLI tool: `qcxms2-to-mgf peaks.dat -o spectrum.mgf`
- [ ] Write unit tests

---

### V2: DreaMS-Based Validation Metric
**Status**: Not Started  
**Impact**: Better similarity than cosine  
**Dependencies**: V1 (MGF converter), DreaMS installation

Use DreaMS embeddings for more meaningful calculated-vs-experimental comparison.

**Concept**:
```python
# benchmark/dreams_metrics.py

from dreams.api import dreams_embeddings
import numpy as np

def dreams_similarity(calculated_mgf: str, experimental_mgf: str) -> float:
    """
    Compare spectra using DreaMS learned embeddings.
    
    Returns cosine similarity in DreaMS embedding space (1024-dim).
    """
    calc_emb = dreams_embeddings(calculated_mgf)
    exp_emb = dreams_embeddings(experimental_mgf)
    
    similarity = np.dot(calc_emb, exp_emb) / (
        np.linalg.norm(calc_emb) * np.linalg.norm(exp_emb)
    )
    return float(similarity)

def dreams_library_search(calculated_mgf: str, library_mgf: str, top_k: int = 5):
    """
    Search for similar experimental spectra using DreaMS.
    """
    from dreams_utils import dreams_search
    return dreams_search(calculated_mgf, library_mgf, top_k=top_k)
```

**Tasks**:
- [ ] Create `benchmark/dreams_metrics.py`
- [ ] Integrate with existing `QCxMS2ValidationMetrics` class
- [ ] Compare DreaMS similarity vs our custom metrics on benzene
- [ ] Benchmark on caffeine, ethanol, toluene
- [ ] Document interpretation of DreaMS similarity scores

---

### V3: Automated Benchmark Pipeline
**Status**: Partial (have components, need integration)  
**Impact**: Systematic validation  
**Dependencies**: V1, database connectors

End-to-end pipeline: molecule -> QCxMS2 -> compare to MassBank/MassSpecGym.

**Workflow**:
```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  Input      │     │   QCxMS2    │     │   Convert   │     │  Compare    │
│  Molecule   │────▶│   Calculate │────▶│   to MGF    │────▶│  to DB      │
│  (SMILES)   │     │   Spectrum  │     │             │     │             │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
                                                                   │
                                                                   ▼
                                                            ┌─────────────┐
                                                            │  Metrics    │
                                                            │  Report     │
                                                            └─────────────┘
```

**Tasks**:
- [ ] Create `benchmark/run_benchmark.py` CLI tool
- [ ] Add molecule-to-XYZ conversion (RDKit/Open Babel)
- [ ] Implement automatic MassBank spectrum fetching
- [ ] Generate HTML/Markdown report with plots
- [ ] Add SLURM job template for batch benchmarking

---

### V4: Multi-Molecule Benchmark Suite
**Status**: Not Started  
**Impact**: Comprehensive validation  
**Dependencies**: V3

Benchmark across diverse molecular classes.

**Proposed Test Set**:
| Molecule | Formula | MW | Class | MassBank Entries |
|----------|---------|-----|-------|------------------|
| Benzene | C6H6 | 78 | Aromatic | 50+ |
| Caffeine | C8H10N4O2 | 194 | Alkaloid | 261 (MassSpecGym) |
| Ethanol | C2H6O | 46 | Alcohol | 30+ |
| Acetone | C3H6O | 58 | Ketone | 25+ |
| Toluene | C7H8 | 92 | Aromatic | 40+ |
| Aspirin | C9H8O4 | 180 | Pharmaceutical | 100+ |
| Glucose | C6H12O6 | 180 | Carbohydrate | 50+ |
| Cholesterol | C27H46O | 386 | Steroid | 20+ |

**Tasks**:
- [ ] Curate test set with XYZ structures
- [ ] Download reference spectra from MassBank
- [ ] Run QCxMS2 calculations (SLURM array job)
- [ ] Compile results into benchmark report
- [ ] Identify systematic errors/missing fragments

---

## Phase 3: ML-Accelerated QCxMS2 (2-4 Weeks)

### A3: Fragment Prioritization via DreaMS
**Status**: Not Started  
**Impact**: 5-10x speedup  
**Dependencies**: DreaMS installation, fragment database

Use DreaMS to predict which fragments are likely, prioritize those in CREST.

**Concept**:
```python
# acceleration/fragment_prioritization.py

from dreams.api import DreaMSSearch

def prioritize_fragments(parent_smiles: str, dreams_atlas: DreaMSSearch) -> list:
    """
    Query DreaMS database for molecules similar to parent.
    Extract common fragment m/z values from neighbors.
    
    Returns list of (m/z, priority_score) tuples.
    """
    # Find similar molecules in 201M spectra database
    similar = dreams_atlas.search_by_structure(parent_smiles, top_k=100)
    
    # Extract common fragment peaks
    fragment_counts = {}
    for spectrum in similar:
        for peak in spectrum.peaks:
            mz_rounded = round(peak.mz, 1)
            fragment_counts[mz_rounded] = fragment_counts.get(mz_rounded, 0) + 1
    
    # Return prioritized by frequency
    return sorted(fragment_counts.items(), key=lambda x: -x[1])
```

**Tasks**:
- [ ] Create `acceleration/fragment_prioritization.py`
- [ ] Implement structure-to-spectrum search in DreaMS
- [ ] Create CREST constraint file generator
- [ ] Benchmark prioritized vs exhaustive search
- [ ] Integrate with QCxMS2 workflow

---

### A4: Formula Pruning via MEDUSA LSTM
**Status**: Not Started  
**Impact**: 3-10x speedup  
**Dependencies**: MEDUSA installation

Use MEDUSA's LSTM to predict molecular formulas, constrain fragmentation search.

**Concept**:
```python
# acceleration/formula_pruning.py

from mass_automation.formula_prediction import predict_formula

def get_formula_constraints(experimental_spectrum: str) -> list:
    """
    Use MEDUSA LSTM to predict formulas for experimental peaks.
    Returns list of formulas to constrain CREST search.
    """
    from mass_automation.experiment import Experiment
    
    exp = Experiment(experimental_spectrum)
    spectrum = exp.summarize()
    
    formulas = []
    for peak in spectrum.peaks:
        predicted = predict_formula(peak.mz, peak.intensity)
        formulas.append(predicted)
    
    return formulas
```

**Tasks**:
- [ ] Create `acceleration/formula_pruning.py`
- [ ] Test MEDUSA LSTM formula prediction accuracy
- [ ] Generate CREST `-formula` constraint files
- [ ] Benchmark constrained vs unconstrained search
- [ ] Handle false negatives (missing formula predictions)

---

### A5: Early Stopping via Embedding Similarity
**Status**: Not Started  
**Impact**: 2-5x speedup  
**Dependencies**: V1, V2

Stop QCxMS2 when partial spectrum is "good enough".

**Concept**:
```python
# acceleration/early_stopping.py

def qcxms2_with_early_stopping(
    input_xyz: str,
    experimental_mgf: str,
    similarity_threshold: float = 0.95,
    max_levels: int = 6
):
    """
    Run QCxMS2 with early stopping based on DreaMS similarity.
    """
    exp_embedding = dreams_embeddings(experimental_mgf)
    
    for level in range(1, max_levels + 1):
        # Run one fragmentation level
        run_qcxms2_level(input_xyz, level)
        
        # Convert partial result to MGF
        partial_mgf = peaks_to_mgf(f"level_{level}/peaks.dat")
        partial_embedding = dreams_embeddings(partial_mgf)
        
        # Check similarity
        similarity = cosine_similarity(partial_embedding, exp_embedding)
        print(f"Level {level}: DreaMS similarity = {similarity:.3f}")
        
        if similarity > similarity_threshold:
            print(f"Early stop at level {level}!")
            break
    
    return partial_mgf
```

**Tasks**:
- [ ] Understand QCxMS2 level-by-level output structure
- [ ] Create `acceleration/early_stopping.py`
- [ ] Implement partial spectrum extraction
- [ ] Benchmark convergence rate on test molecules
- [ ] Determine optimal similarity threshold

---

### A6: Pre-computed Fragment Database
**Status**: Not Started  
**Impact**: 10-100x speedup (long-term)  
**Dependencies**: A2 (cache system)

Build shared database of fragments with pre-computed NEB barriers.

**Schema**:
```sql
CREATE TABLE fragments (
    id INTEGER PRIMARY KEY,
    mz REAL NOT NULL,
    formula TEXT NOT NULL,
    smiles TEXT NOT NULL,
    inchi_key TEXT NOT NULL,
    barrier_ev REAL,
    parent_smiles TEXT,
    calculation_level TEXT,
    source TEXT,  -- 'qcxms2', 'literature', 'ml_predicted'
    created_at TIMESTAMP
);

CREATE INDEX idx_mz ON fragments(mz);
CREATE INDEX idx_inchi ON fragments(inchi_key);
```

**Tasks**:
- [ ] Design database schema
- [ ] Create `database_integration/fragment_db/` module
- [ ] Implement fragment fingerprint matching
- [ ] Seed with fragments from completed calculations
- [ ] Add barrier lookup to QCxMS2 NEB step
- [ ] Create web API for fragment queries (optional)

---

## Phase 4: Research Directions (Future)

### R1: ML-Predicted NEB Barriers
**Status**: Research  
**Impact**: 100-1000x speedup  
**Dependencies**: Training data from QCxMS2 calculations

Train neural network to predict transition state barriers.

**Approach**:
```
Training Data:
- Input: (Reactant SMILES, Product SMILES, Bond broken)
- Output: ΔG‡ (barrier height in eV)
- Source: All NEB calculations from QCxMS2 runs

Model Options:
- Graph Neural Network (GNN) on molecular graphs
- Transformer on SMILES strings
- SchNet/PaiNN on 3D structures

Validation:
- Hold-out test set from known barriers
- Compare to DFT/GFN2 reference values
```

**Tasks**:
- [ ] Collect NEB barrier data from existing calculations
- [ ] Design training dataset format
- [ ] Implement GNN baseline model
- [ ] Train and validate on test set
- [ ] Integrate with QCxMS2 as optional fast mode

---

### R2: GNN Direct Fragmentation Prediction
**Status**: Research  
**Impact**: 1000x speedup (screening)  
**Dependencies**: Large training dataset

Skip QCxMS2 entirely - predict spectrum directly from structure.

**Approach**:
```
Training Data:
- Input: Molecular graph + DreaMS embeddings
- Output: (m/z, intensity) peaks
- Source: MassBank + MassSpecGym + QCxMS2 calculations

Model:
- Graph attention network for fragmentation sites
- Set Transformer for peak prediction
- Use DreaMS pre-trained embeddings as features

Use Case:
- Fast screening of large compound libraries
- Refine top candidates with full QCxMS2
```

**Tasks**:
- [ ] Literature review of existing ML fragmentation models
- [ ] Collect training data from public databases
- [ ] Implement baseline model
- [ ] Compare accuracy vs QCxMS2
- [ ] Publish if results are promising

---

## Completed Work

### Fragment Cache System (January 15, 2026)
- [x] SQLite-based fragment cache (`cache/fragment_cache.py`)
- [x] QCxMS2 output parser (`cache/qcxms2_output_parser.py`)
- [x] CLI tool for cache management (`cache/cli.py`)
- [x] Production database seeded with 132 fragments from benzene calculations
- [x] Database location: `/mnt/beegfs/software/qcxms2/fragment_cache.db`

**Usage**:
```python
from database_integration.cache import FragmentCache

cache = FragmentCache()  # Uses default database
result = cache.find_fragment(mz=77.104, formula="C6H5")
if result:
    print(f"Cached barrier: {result.barrier_kcal_mol} kcal/mol")
```

### Database Integration Module (January 2026)
- [x] MassBank v3 API connector (`connectors/massbank.py`)
- [x] MassSpecGym connector (`connectors/massspecgym.py`)
- [x] QCxMS2 parser (`parsers/qcxms2_parser.py`)
- [x] Benchmark suite (`benchmark/benchmark_suite.py`)
- [x] QCxMS2-specific validation metrics (`benchmark/metrics.py`)
  - [x] `calculate_fragmentation_pattern_match()`
  - [x] `calculate_mass_accuracy_score()`
  - [x] `calculate_intensity_rank_correlation()`
  - [x] `calculate_qcxms2_validation_metrics()`

### IEE Parameter Optimization (January 2026)
- [x] IEE sweep job script (`run_iee_mlip_sweep.slurm`)
- [x] Tested ieeatm values: 0.4, 0.6, 0.8, 1.0
- [x] Found optimal ieeatm=0.4 for benzene
- [x] Documented in `VALIDATION_REPORT.md`

### MLIP Integration (January 2026)
- [x] MLIP-accelerated NEB calculations working
- [x] AIMNet2 server at id-gpu01:8888
- [x] CREST MLIP wrapper script

---

## Known Issues

### Missing Fragmentation Pathways
**Problem**: QCxMS2 does not find m/z 77, 51, 39 for benzene even though these are major experimental peaks.

**Possible Causes**:
1. CREST msreact not exploring these reaction channels
2. NEB failing to find transition states
3. Barriers too high relative to IEE

**Investigation Needed**:
- [ ] Check CREST msreact output for phenyl cation (m/z 77)
- [ ] Try different `-msnbonds` and `-msnshifts` settings
- [ ] Use higher-level TS optimization (`-tslevel r2scan3c`)

### Intensity Prediction Accuracy
**Problem**: Relative intensities differ significantly from experiment.

**Current Understanding**:
- Eyring/RRKM theory is correct approach
- IEE parameter strongly affects results
- Barrier heights from GFN2/MLIP may be inaccurate

**Mitigation**:
- Lower ieeatm (0.4) gives better results
- Consider ML barrier prediction (R1)

### Caffeine CID Timeout
**Problem**: Caffeine CID calculation times out at 4 hours.

**Solutions**:
- [ ] Run with 12-hour time limit
- [ ] Use MLIP acceleration (`-mlip`)
- [ ] Reduce fragmentation depth (`-nfrag 3`)

---

## File Index

| Path | Description |
|------|-------------|
| `database_integration/` | Database connectors and validation tools |
| `database_integration/connectors/` | MassBank, MassSpecGym connectors |
| `database_integration/parsers/` | QCxMS2 output parsers |
| `database_integration/benchmark/` | Benchmark suite and metrics |
| `database_integration/converters/` | Format converters (planned) |
| `database_integration/cache/` | Fragment cache (planned) |
| `acceleration/` | ML acceleration modules (planned) |

---

## References

### QCxMS2
- Paper: Gorges & Grimme, PCCP 27, 6899-6911 (2025)
- DOI: [10.1039/D5CP00316D](https://doi.org/10.1039/D5CP00316D)

### DreaMS
- Paper: Bushuiev et al., Nature Biotechnology (2025)
- DOI: [10.1038/s41587-025-02663-3](https://doi.org/10.1038/s41587-025-02663-3)
- Cluster setup: `/mnt/beegfs/software/dreams/setup-dreams.sh`

### MEDUSA
- Paper: Abate-Pella et al. (2024)
- Cluster setup: `/mnt/beegfs/software/medusa/setup-medusa.sh`
- Documentation: `/mnt/beegfs/software/medusa/AI_AGENT_COMPLETE_GUIDE.md`

### MassSpecGym
- Paper: NeurIPS 2024 Datasets and Benchmarks
- Dataset: 231K MS/MS spectra

### MLIP Servers
- AIMNet2 JSON API: `id-gpu01.materials.local.cmu.edu:8888`
- MLIP-G16 API: `gpg-boltzmann:5003`
