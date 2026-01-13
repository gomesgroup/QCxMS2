# Mass Spectrometry Databases Research

**Date**: January 13, 2026  
**Purpose**: Identify additional MS databases for QCxMS2 integration

---

## Currently Implemented

| Database | Status | Records | API | Notes |
|----------|--------|---------|-----|-------|
| **MassBank EU** | ✅ Full | 40K+ | REST | Free, open access |
| **MoNA** | ⚠️ Partial | 500K+ | REST | Free, open access |
| **NIST MS** | ⚠️ Stub | Large | None | Commercial, no public API |
| **SDBS** | ⚠️ Stub | 34K | None | 50/day limit |

---

## High Priority - Should Implement

### 1. mzCloud (HighChem)
**URL**: https://www.mzcloud.org  
**Records**: 32,330 compounds, 51,002 trees, 16.5M+ spectra  
**API**: REST API available  
**Access**: Free with login (Google/Facebook)  
**Data Types**: MSn trees, high-resolution MS/MS  
**Ionization**: ESI+, ESI-  
**Key Features**:
- Spectral trees (MSn data)
- Structure search
- Substructure search
- m/z search
- Tree search for confident identification

**Integration Priority**: ⭐⭐⭐⭐⭐ HIGH
**Reason**: Large dataset, MSn trees, REST API

---

### 2. GNPS (Global Natural Products Social Molecular Networking)
**URL**: https://gnps.ucsd.edu  
**Records**: Massive (natural products focus)  
**API**: REST API available  
**Access**: Free, open access  
**Data Types**: MS/MS molecular networking  
**Key Features**:
- Molecular networking
- Spectral library search
- Community-contributed spectra
- Integration with MassIVE repository

**Integration Priority**: ⭐⭐⭐⭐⭐ HIGH
**Reason**: Large community database, open access, REST API

---

### 3. LIPID MAPS
**URL**: https://www.lipidmaps.org  
**Records**: Comprehensive lipid spectra  
**API**: REST API + SPARQL endpoint  
**Access**: Free, open access (CC BY 4.0)  
**Data Types**: Lipid MS/MS spectra  
**Key Features**:
- Lipid Structure Database (LMSD)
- Standard Spectra Database (1,800+ spectra)
- MS Analysis tools
- Bulk search capability
- Oxylipin and sterol standards

**Integration Priority**: ⭐⭐⭐⭐ HIGH
**Reason**: Specialized lipid database, REST API, well-documented

---

### 4. Metabolomics Workbench / NMDR
**URL**: https://www.metabolomicsworkbench.org  
**Records**: 4,441 studies, 4,000+ public  
**API**: REST API available  
**Access**: Free, NIH-sponsored  
**Data Types**: LC-MS, GC-MS, NMR  
**Key Features**:
- National Metabolomics Data Repository
- RefMet standardized nomenclature
- 4.5M+ m/z features searchable
- Raw data availability
- Study-level metadata

**Integration Priority**: ⭐⭐⭐⭐ HIGH
**Reason**: NIH resource, comprehensive, REST API

---

### 5. MetaboLights (EMBL-EBI)
**URL**: https://www.ebi.ac.uk/metabolights/  
**Records**: Large repository  
**API**: REST API available  
**Access**: Free, open access  
**Data Types**: Metabolomics experiments  
**Key Features**:
- EBI/EMBL-hosted
- Cross-species data
- Reference spectra
- Endorsed by major journals
- Compound library

**Integration Priority**: ⭐⭐⭐⭐ HIGH
**Reason**: Major repository, REST API, journal endorsement

---

## Medium Priority

### 6. CFM-ID (Competitive Fragmentation Modeling)
**URL**: https://cfmid.wishartlab.com  
**Records**: Predicted spectra (not experimental)  
**API**: Web API  
**Access**: Free  
**Data Types**: Predicted MS/MS spectra  
**Key Features**:
- Spectra prediction from structure
- Peak assignment
- Compound identification
- EI-MS and ESI-MS/MS prediction

**Integration Priority**: ⭐⭐⭐ MEDIUM
**Reason**: Predicted (not experimental), but useful for validation

---

### 7. HMDB (Human Metabolome Database)
**URL**: https://hmdb.ca  
**Records**: 220,000+ metabolites  
**API**: REST API available  
**Access**: Free  
**Data Types**: MS/MS, NMR  
**Key Features**:
- Human metabolome focus
- Extensive metadata
- Biofluid concentrations
- Disease associations

**Integration Priority**: ⭐⭐⭐ MEDIUM
**Reason**: Human-focused, may have limited relevance for general MS

---

### 8. METLIN (Scripps)
**URL**: https://metlin.scripps.edu  
**Records**: 1M+ metabolites  
**API**: Limited (requires subscription for full access)  
**Access**: Free basic, premium for full  
**Data Types**: MS/MS  
**Key Features**:
- Large metabolite database
- In-silico fragmentation
- Isotope patterns

**Integration Priority**: ⭐⭐⭐ MEDIUM
**Reason**: Limited free API access

---

### 9. SpectraBase (Wiley)
**URL**: https://spectrabase.com  
**Records**: Millions of spectra  
**API**: Limited (10 searches/30 days free)  
**Access**: Free limited, subscription for full  
**Data Types**: MS, NMR, IR, Raman, UV-Vis  
**Key Features**:
- Multi-spectral database
- Commercial-grade data

**Integration Priority**: ⭐⭐ LOW-MEDIUM
**Reason**: Strict usage limits, commercial focus

---

## Lower Priority / Specialized

### 10. ChEBI (Chemical Entities of Biological Interest)
**URL**: https://www.ebi.ac.uk/chebi/  
**API**: REST + SPARQL  
**Focus**: Chemical ontology, some MS data  
**Integration Priority**: ⭐⭐ LOW

### 11. PubChem (already in Explorer)
**URL**: https://pubchem.ncbi.nlm.nih.gov  
**API**: REST API  
**Focus**: General chemical data, some MS  
**Integration Priority**: ⭐⭐ LOW (already available via Explorer)

### 12. ReSpect (RIKEN Plant MS/MS)
**URL**: http://spectra.psc.riken.jp/  
**Focus**: Plant metabolites  
**Integration Priority**: ⭐⭐ LOW

### 13. MassIVE (UCSD)
**URL**: https://massive.ucsd.edu  
**Focus**: Raw MS data repository  
**Integration Priority**: ⭐⭐ LOW (raw data, not library)

---

## Summary: Recommended Implementation Order

### Phase 1 (Immediate)
1. **mzCloud** - Large dataset, MSn trees, good API
2. **GNPS** - Community database, molecular networking
3. **LIPID MAPS** - Specialized lipids, REST API

### Phase 2 (Short-term)
4. **Metabolomics Workbench** - NIH resource, comprehensive
5. **MetaboLights** - EBI repository, journal standard

### Phase 3 (Medium-term)
6. **CFM-ID** - Predicted spectra for validation
7. **HMDB** - Human metabolome
8. **METLIN** - If API access improves

---

## API Comparison

| Database | API Type | Auth Required | Rate Limit | Documentation |
|----------|----------|---------------|------------|---------------|
| MassBank | REST | No | ~2/s | Good |
| MoNA | REST | No | ~2/s | Good |
| mzCloud | REST | Yes (free) | Unknown | Limited |
| GNPS | REST | No | Unknown | Good |
| LIPID MAPS | REST+SPARQL | No | Unknown | Excellent |
| MetaboLights | REST | No | Unknown | Good |
| Metabolomics WB | REST | No | Unknown | Good |
| HMDB | REST | No | Unknown | Good |
| CFM-ID | Web | No | Unknown | Good |

---

## Data Format Considerations

Most databases use one of these formats:
- **MSP** (NIST format) - Most common
- **MGF** (Mascot Generic Format) - Common for MS/MS
- **JSON** - Modern REST APIs
- **mzML/mzXML** - Raw data formats (not library)

Our `utils.py` already supports MSP, MGF, and TXT parsing.

---

## EI-MS Specific Databases

For QCxMS2's EI-MS focus, the most relevant databases are:
1. **MassBank** (has EI spectra) ✅
2. **NIST MS** (gold standard for EI, but no API)
3. **mzCloud** (has some EI data)
4. **Metabolomics Workbench** (has GC-MS/EI data)

---

## Next Steps

1. Implement **mzCloud** connector (highest value addition)
2. Implement **GNPS** connector (large community database)
3. Implement **LIPID MAPS** connector (specialized, good API)
4. Test and validate existing connectors
5. Create benchmark comparison suite

---

## References

1. mzCloud: https://www.mzcloud.org/About
2. GNPS: Wang et al., Nat Biotechnol. 2016;34(8):828-837
3. LIPID MAPS: Sud et al., Nucleic Acids Res. 2007;35:D527-D532
4. MetaboLights: Haug et al., Nucleic Acids Res. 2020;48(D1):D440-D444
5. Metabolomics Workbench: Sud et al., Nucleic Acids Res. 2016;44(D1):D463-D470
6. CFM-ID: Wang et al., Nucleic Acids Res. 2022;50(W1):W165-W174
