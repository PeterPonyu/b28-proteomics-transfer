# Third-party boundaries

This export licenses only Zeyu Fu's code and authored research content.

- **CPTAC and UMICH data:** processed reference-normalized protein abundance
  tables, clinical/sample mapping inputs and source annotations are third-party
  data. They are excluded; filenames/checksums and public scientific identifiers
  are retained as provenance. The data source's access, use and citation terms
  govern acquisition and reuse. `DATA_SOURCES.md` distinguishes the actual
  recorded source from unrecorded PDC accession/version fields.
- **CPTAC Python package:** a dependency (`cptac 1.5.14`), not vendored here. Its
  source carries Apache-2.0 notices. Its upstream data repository identifier is
  not this B28 research-release DOI. No CPTAC downloader or credential is copied.
- **Other Python/R/TeX dependencies:** NumPy, pandas, SciPy, scikit-learn,
  XGBoost, Matplotlib, PyMuPDF, Pillow, ggplot2 and the other recorded R packages,
  Poppler, TeX Live and the ACS `achemso` document class are externally supplied
  dependencies. Their own licenses govern them; this repository does not
  redistribute them or claim MIT/CC BY ownership of them.
- **Fonts:** current published figures embed or outline rendering glyphs; exact
  reconstruction expects separately licensed Arial regular/bold. No standalone
  Arial font files are included. TeX fonts remain part of the external TeX
  installation. A missing required font fails rather than silently substituting
  and asserting equivalent output.
- **Bibliography:** bibliographic facts, titles and DOIs attribute the cited
  authors and publishers; this release does not redistribute cited full-text
  papers or relicense them.

Historical manifest references to private submission paths do not include those
files. A retained historical OV summary is author-generated contextual output,
not the original OV abundance or annotation file and not a new prediction replay.
