# Actual B28 data sources and reuse boundary

The six current cohorts were read through `cptac 1.5.14`, source `umich`,
using `Ccrcc`, `Luad`, `Ucec`, `Hnscc`, `Lscc`, and `Pdac` objects.
The exact frozen fit-time abundance hashes come from the included
`analysis/revision_20260922/candidate-v1/candidate-results.json`.
Source values are UMICH reference-intensity-subtracted abundance, not
raw mass-spectrometer total signal. The stable ENSP/gene mapping is
published separately in the corrected accession audit.

The existing installed CPTAC downloader records upstream **data DOI
10.5281/zenodo.8394329**, index DOI 10.5281/zenodo.7897498 and an index
download record 10573663. These are third-party source identifiers, not
B28 release identifiers. They were inspected locally, not reconstructed
from a new download. The cached package/index and frozen hashes, rather
than a guessed release label, define the exact observed input.

**PDC provenance limit:** CPTAC data may be accessible through its
Proteomic Data Commons ecosystem, but this fit actually used the UMICH
tables exposed by the CPTAC package/Zenodo archive. The existing fit
records contain no PDC study accession or separately named per-cohort
data version. Those fields are explicitly null in `provenance/data-sources.json`.
No PDC accession is invented or assigned from a cancer-type label.

## Frozen abundance files

All six local basenames are `Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
The cohort directory and prefixed upstream key distinguish them. The following
are provenance descriptions, not redistributed third-party files.

### Ccrcc - Clear cell renal cell carcinoma

- Local relative file: `umich-ccrcc/Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Upstream repository file key: `umich-ccrcc-proteomics-Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Frozen fit-time SHA-256: `fd111b453c205ce2f006d4d0ea9403b5e14e492df88e6611d1694c1e583fb6f4`; 19,335,367 bytes.
- CPTAC package: 1.5.14; source UMICH; PDC accession and separate cohort version not recorded.
- Current cached-index MD5 matches frozen local gzip: `True`.

### Luad - Lung adenocarcinoma

- Local relative file: `umich-luad/Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Upstream repository file key: `umich-luad-proteomics-Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Frozen fit-time SHA-256: `c14afab469226dc7fb709109cbcc0228a9189b9cd6660519573714cf0b3f3bba`; 23,998,218 bytes.
- CPTAC package: 1.5.14; source UMICH; PDC accession and separate cohort version not recorded.
- Current cached-index MD5 matches frozen local gzip: `True`.

### Ucec - Uterine corpus endometrial carcinoma

- Local relative file: `umich-ucec/Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Upstream repository file key: `umich-ucec-proteomics-Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Frozen fit-time SHA-256: `c050e457a521bbbf36b28a0f556b4a530403eea98317673b78e94f676cee91cf`; 15,733,938 bytes.
- CPTAC package: 1.5.14; source UMICH; PDC accession and separate cohort version not recorded.
- Current cached-index MD5 matches frozen local gzip: `True`.

### Hnscc - Head and neck squamous cell carcinoma

- Local relative file: `umich-hnscc/Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Upstream repository file key: `umich-hnscc-proteomics-Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Frozen fit-time SHA-256: `b2f2396ae04c569282c1792625a433b517854dff42445c3655257291455dc851`; 19,369,270 bytes.
- CPTAC package: 1.5.14; source UMICH; PDC accession and separate cohort version not recorded.
- Current cached-index MD5 matches frozen local gzip: `True`.

### Lscc - Lung squamous cell carcinoma

- Local relative file: `umich-lscc/Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Upstream repository file key: `umich-lscc-proteomics-Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Frozen fit-time SHA-256: `24370f80b91ece04a5a4d8b6099fcf1155a0f24ef76f2614890aa2a15b8b99b0`; 24,449,704 bytes.
- CPTAC package: 1.5.14; source UMICH; PDC accession and separate cohort version not recorded.
- Current cached-index MD5 matches frozen local gzip: `True`.

### Pdac - Pancreatic ductal adenocarcinoma

- Local relative file: `umich-pdac/Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Upstream repository file key: `umich-pdac-proteomics-Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`.
- Frozen fit-time SHA-256: `336d148faa5ca48054acdb6c03c723b6da9873eaeb6b641b028390ecaaad29a5`; 20,504,572 bytes.
- CPTAC package: 1.5.14; source UMICH; PDC accession and separate cohort version not recorded.
- Current cached-index MD5 matches frozen local gzip: `True`.

## Loader dependencies and exact-version caution

For CCRCC, LUAD, UCEC, LSCC and PDAC the UMICH loaders additionally consume
`aliquot_to_patient_ID.tsv.gz`. Their cohort-prefixed repository keys and
current observed SHA-256 values are recorded in `provenance/data-sources.json`.
The original fit did not hash these mappings independently; their present
observed hashes are **not** retroactive evidence of fit-time identity.
HNSCC instead uses source-column IDs, duplicate handling and `.N`/`.C`
suffix parsing in its loader. The public corrected sample/fold maps are
the direct stored-output identity reference.

The cached index checksums do not necessarily match the frozen local gzip
bytes. Do not silently redownload newer files and call them the reviewed
matrix. Match the recorded SHA-256, source-file inventory, sample labels,
patient identity and protein mapping before fitting. This export neither
claims a new upstream download nor redistributes abundance/mapping inputs.

## Historical OV context

SI S7 refers only to the narrowly retained historical summary
`provenance/historical/ov-summary-20260911.json`. Its independent UMICH OV
abundance file and `OV_sample_TMT_annotation_UMich_GENCODE34_0315.csv.gz`
have original source hashes/bytes in that record but are omitted here.
The prior source/annotation check is historical. No prediction-level OV
replay, new score or independence/deconfounding claim is made.

## Acquisition, attribution and licensing

Obtain upstream files from the original CPTAC/UMICH distribution, respecting
its access/use/citation terms. The CPTAC package documentation is
https://paynelab.github.io/cptac/ and its package paper DOI is
10.1021/acs.jproteome.0c00919 (already cited in the included bibliography).
The current manuscript also cites the original CPTAC/pan-cancer sources.
The B28 MIT/CC BY licenses do not replace upstream data terms.

Public sample/patient/protein identifiers in derived maps preserve traceability,
not private clinical records or raw quantitative proteomics. Full model fitting
requires separately acquired upstream matrices; stored prediction/metric and
offline PDF/figure verification do not. No fit or download was performed
during this export.
