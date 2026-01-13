# QCxMS2 Database Integration

Modular system for accessing experimental mass spectrometry databases to compare with calculated QCxMS2 spectra.

**Architecture inspired by Explorer's database_integration system.**

## Features

- **Modular Design**: Easy to add new database connectors
- **Multiple Databases**: MassBank, MoNA, NIST MS (stub), SDBS (stub)
- **REST API Access**: Programmatic access to 40K+ experimental spectra
- **Spectrum Comparison**: Cosine similarity, dot product, Euclidean distance
- **Caching**: Automatic result caching to minimize API calls
- **Rate Limiting**: Respects database usage policies

## Supported Databases

| Database | Status | API | Records | Notes |
|----------|--------|-----|---------|-------|
| **MassBank EU** | ✅ Full | REST | 40K+ | Free, no API key |
| **MoNA** | ✅ Partial | REST | 500K+ | Free, no API key |
| **NIST MS** | ⚠️ Stub | None | Large | No public API, web scraping needed |
| **SDBS** | ⚠️ Stub | None | ~34K | 50 queries/day limit, no official API |

## Installation

### Dependencies

```bash
# Install Python dependencies
pip install requests numpy

# Or use QCxMS2 environment
source /mnt/beegfs/software/qcxms2/setup-qcxms2.sh
```

### Add to Python Path

```python
import sys
sys.path.insert(0, "/mnt/beegfs/software/qcxms2")

from database_integration import create_ms_connector, DatabaseType
```

## Quick Start

### Example 1: Search by Molecular Formula

```python
from database_integration import (
    create_ms_connector,
    DatabaseType,
    SpectrumSearchCriteria,
    IonizationMode,
)

# Create MassBank connector
massbank = create_ms_connector(DatabaseType.MASSBANK)

# Search for glucose EI spectra
criteria = SpectrumSearchCriteria(
    formula="C6H12O6",
    ionization_mode=IonizationMode.EI,
    max_results=10,
)

results = massbank.search_spectrum(criteria)

if results.success:
    print(f"Found {results.total_results} spectra")
    for spectrum in results.spectra:
        print(f"  {spectrum.spectrum_id}: {spectrum.name}")
        print(f"    Peaks: {len(spectrum.peaks)}")
        print(f"    URL: {spectrum.url}")
```

### Example 2: Compare Calculated vs. Experimental

```python
from database_integration import compare_spectra, SimilarityMetric, Peak

# Your calculated QCxMS2 peaks
calculated_peaks = [
    Peak(mz=43.0, intensity=100.0),
    Peak(mz=57.0, intensity=80.0),
    Peak(mz=71.0, intensity=60.0),
]

# Experimental peaks from database
experimental_spectrum = massbank.get_spectrum_by_id("MSBNK-BAFG-CSL2328300")

# Compare spectra
comparison = compare_spectra(
    calculated_peaks,
    experimental_spectrum.peaks,
    metric=SimilarityMetric.COSINE,
    mz_tolerance=0.5,
)

print(f"Similarity: {comparison.similarity_score:.3f}")
print(f"Matched peaks: {comparison.matched_peaks}/{comparison.total_peaks_query}")
```

### Example 3: Batch Comparison Against Database

```python
from database_integration import BatchComparisonResult

# Search for all glucose spectra
criteria = SpectrumSearchCriteria(formula="C6H12O6", max_results=50)
results = massbank.search_spectrum(criteria)

# Compare your calculated spectrum against all results
comparisons = []
for exp_spectrum in results.spectra:
    comp = compare_spectra(
        calculated_peaks,
        exp_spectrum.peaks,
        query_id="calculated",
        reference_id=exp_spectrum.spectrum_id,
    )
    comparisons.append(comp)

# Create batch result
batch = BatchComparisonResult(
    query_spectrum_id="calculated",
    comparisons=comparisons,
)

# Get top 5 matches
top_matches = batch.get_top_n(5)
for i, match in enumerate(top_matches, 1):
    print(f"{i}. {match.reference_spectrum_id}: {match.similarity_score:.3f}")
```

## Architecture

The system follows Explorer's modular database integration design:

```
database_integration/
├── __init__.py          # Public API
├── enums.py             # DatabaseType, IonizationMode, SimilarityMetric
├── data_classes.py      # Peak, ExperimentalSpectrum, MSSearchResult
├── search_criteria.py   # SpectrumSearchCriteria
├── config.py            # DatabaseConfig with defaults
├── base.py              # MSDatabaseConnector (abstract base)
├── factory.py           # create_ms_connector()
├── utils.py             # compare_spectra(), parse_spectrum_file()
└── connectors/
    ├── __init__.py
    ├── massbank.py      # MassBankConnector (full implementation)
    ├── mona.py          # MoNAConnector (partial)
    ├── nist_ms.py       # NISTMSConnector (stub)
    └── sdbs.py          # SDBSConnector (stub)
```

### Key Design Patterns

1. **Factory Pattern**: `create_ms_connector(DatabaseType.MASSBANK)`
2. **Strategy Pattern**: Different similarity metrics (cosine, dot product, etc.)
3. **Template Method**: `MSDatabaseConnector` base class with common HTTP/caching logic
4. **Data Classes**: Type-safe data structures with validation

## Integration with QCxMS2

### Option 1: Post-Processing Script

```python
#!/usr/bin/env python3
"""
Compare QCxMS2 output with experimental database.
"""
import sys
sys.path.insert(0, "/mnt/beegfs/software/qcxms2")

from database_integration import *

# Parse QCxMS2 output
qcxms_peaks = parse_qcxms_output("qcxms.out")

# Search database
massbank = create_ms_connector(DatabaseType.MASSBANK)
results = massbank.search_spectrum(SpectrumSearchCriteria(formula="C10H16"))

# Compare and rank
for exp_spec in results.spectra:
    sim = compare_spectra(qcxms_peaks, exp_spec.peaks)
    print(f"{exp_spec.name}: {sim.similarity_score:.3f}")
```

### Option 2: Direct Integration (Future)

```fortran
! In QCxMS2 Fortran code, call Python via subprocess or f2py
! This would require modifications to QCxMS2 source
```

## API Usage Limits

| Database | Limit | Notes |
|----------|-------|-------|
| MassBank | ~2 req/s | Conservative, no official limit |
| MoNA | ~2 req/s | Conservative, no official limit |
| NIST MS | N/A | No API (web scraping) |
| SDBS | 50/day | Strict limit, use sparingly |

The system automatically respects rate limits via `DatabaseConfig.rate_limit_per_second`.

## Caching

Results are automatically cached in `~/.cache/qcxms2/ms_databases/` for 24 hours.

To disable caching:

```python
config = DatabaseConfig.for_massbank(cache_enabled=False)
connector = MassBankConnector(config)
```

To change cache TTL:

```python
config = DatabaseConfig.for_massbank(cache_ttl=3600)  # 1 hour
```

## Error Handling

```python
try:
    results = massbank.search_spectrum(criteria)
    if not results.success:
        print(f"Error: {results.error_message}")
except Exception as e:
    print(f"Exception: {e}")
```

## Adding New Databases

To add a new database connector:

1. Create `connectors/mydatabase.py`
2. Inherit from `MSDatabaseConnector`
3. Implement abstract methods:
   - `test_connection()`
   - `search_spectrum(criteria)`
   - `get_spectrum_by_id(id)`
4. Add to `enums.DatabaseType`
5. Add to `factory.create_ms_connector()`

See `connectors/massbank.py` for a complete example.

## Citation

If you use this database integration system, please cite:

**MassBank**:
> Horai et al. MassBank: a public repository for sharing mass spectral data for life sciences. J. Mass Spectrom. 2010, 45(7):703-714.

**MoNA**:
> MassBank of North America (MoNA). https://mona.fiehnlab.ucdavis.edu/

## License

This module is part of QCxMS2 and follows the same license.

## References

- MassBank API: https://github.com/MassBank/MassBank-web
- MoNA API: https://mona.fiehnlab.ucdavis.edu/rest/api-doc.html
- NIST Chemistry WebBook: https://webbook.nist.gov/chemistry/
- SDBS: https://sdbs.db.aist.go.jp/
