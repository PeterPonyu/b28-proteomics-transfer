# B28: accession-aligned proteomics evaluation-protocol contrast

Public research repository: https://github.com/PeterPonyu/b28-proteomics-transfer

This is an allowlisted public export of the current corrected B28 research
analysis, not a journal-submission package or a clinical-validation release.
Scientific status: **SCIENTIFIC_HOLD**. Tests and hashes establish reproducible
file identity and consistency, not scientific clearance, causality, clinical
utility, or independent replication.

Current version DOI: **https://doi.org/10.5281/zenodo.23101188**.
Research release tag: **`author-review-2026-10-02`**, version
`2026.10.02-author-review`. The matching code, current Article/SI, stored
predictions and integrity manifest are distributed together. The original
B28 Zenodo concept is preserved; earlier archives are historical, not this
corrected release. The identifier alone is not proof of publication: remote
publication and anonymous-download verification are recorded separately in the
dated release receipts. The manuscript DAS and current compiled PDFs refer to
this version. No journal submission or scientific clearance is claimed.

## Included current evidence

- The current 24-page Article (`manuscript/manuscript.pdf`) and 6-page supporting
  information (`manuscript/SI.pdf`), their consumed TeX sources, bibliography,
  tables, figures F0--F7, TOC, current graphical abstract, and build scripts.
- The actual corrected analysis programs under `analysis/revision_20260922/`:
  accession/sample/patient reanalysis, matched LOGO sensitivities, recalibration,
  and fixed-prediction conditional patient resampling.
- Stable ENSP accession and gene-mapping audits; all 1,191 sample identities;
  cohort-scoped patient IDs; folds for sample CV, patient CV, and cohort LOGO;
  and the identity-keyed corrected five-fold replay table.
- 7,146 primary/comparator prediction rows, 2,382 matched sensitivity prediction
  rows, 1,191 stored recalibration rows, recorded metrics, and independent
  numerical-verification dependencies. Dated `candidate-v1` and
  `candidate-v2-side-studies` names are preserved for traceability; these are the
  corrected results consumed by the current Article, not obsolete manuscripts.
- A narrowly retained OV summary in `provenance/historical/`, explicitly not
  prediction-replayed or promoted as independently verified current evidence.

The corrected population is 1,191 samples from 719 patients, with 689 tumors,
496 adjacent-normal samples and six HNSCC cored-normal samples; the primary
intersection contains 8,906 stable ENSP accessions. Study/cancer confounding and
target-cohort feature-availability conditioning remain. No new scores, model
fits, inference, uncertainty claims or readiness claims were added during export.

## Reproduction: three distinct boundaries

Run all commands from the repository root with an existing compatible Python
environment. No command below installs dependencies.

1. **Stored predictions, metrics and identity checks (no upstream matrix):**

   ```bash
   python reproducibility/verify_method_contract.py
   python analysis/revision_20260922/verify_candidate.py
   python analysis/revision_20260922/verify_postprocess.py
   python analysis/revision_20260922/verify_side_studies.py
   python -m pytest -q -p no:cacheprovider
   ```

   Default checks do not import CPTAC or read the author's recorded absolute
   raw-data paths. They explicitly report that the upstream matrix was not
   checked. New verification reports go to `output/verification/`; frozen
   scientific predictions and QA files are not overwritten.

2. **Offline figure/manuscript verification (no fit, no PDF edit):**

   ```bash
   bash manuscript/scripts/rebuild_current.sh
   python manuscript/scripts/verify_integration.py
   ```

   The first command checks the public file lock, verifies stored metrics,
   renders all supplied PDF/figure pages using Poppler, checks text and
   geometry, verifies the accepted current compilation inputs, and reuses the prior
   301-cell numerical record only after checking its 39 included dependencies.
   The second independently recomputes those 301 existing metric cells from
   stored rows; it does not predict, resample or refit a model. Rendering output
   is local and disposable. The successful default does not claim a new TeX
   compilation or revalidation of third-party raw data.

   For an **explicit local build**, with the existing Python/R/TeX/font
   requirements described in `environment/README.md`, run:

   ```bash
   bash manuscript/scripts/rebuild_current.sh --build
   ```

   This regenerates current tables/figures and compiles Article/SI. It does not
   change the scientific evidence or refresh the public lock. If bytes differ,
   the final fail-closed lock check reports the changed inputs. A rebuild is not
   automatic permission to bless changed sources or replace a published PDF.

3. **Full model fitting (upstream matrix required, not run for this export):**

   Obtain the exact UMICH abundance/mapping inputs in `DATA_SOURCES.md` and use
   the recorded fit environment, not the offline-verifier environment. The
   programs refuse to overwrite their existing output. A separate output root
   can be selected explicitly:

   ```bash
   export B28_FIT_ROOT=/absolute/path/to/a/new-b28-fit
   python analysis/revision_20260922/reanalyse.py
   python analysis/revision_20260922/side_studies.py
   python analysis/revision_20260922/postprocess.py
   python analysis/revision_20260922/fixed_prediction_sensitivity.py
   ```

   CPTAC may initialize a remote index on import; fitting is not part of the
   network-free verification path. Existing local source tables must be
   populated first. Repeated fits are not promised bit-for-bit identical: the
   original XGBoost settings did not explicitly set a random state. Validate
   new outputs separately and do not relabel them as the stored release.
   `--upstream-root /path/to/cptac/data` on the two stored-metric validators
   optionally checks exact upstream abundance hashes; it fails on absent or
   changed files and does not download them.

## Integrity, history and licensing

`PUBLIC-MANIFEST.json`, `reproducibility/method-contract.json`, and
`SHA256SUMS.txt` describe the actual included public files. Check the independent
transport inventory with `sha256sum -c SHA256SUMS.txt`. Unexpected inputs,
missing files, path traversal, symlinks, changed sources and changed predictions
are rejected. Local generated output, bytecode, Git metadata and TeX scratch
files are not release payload and are excluded from file-lock comparison.

The original capsule manifests are preserved only in `provenance/historical/`.
They reference excluded backups, old rendered pages and unsynchronized private
submission packages, so their dependency closure is **not** the current public
manifest. `provenance/export-source-ledger.json` records original source hashes;
`EXPORT_NOTES.md` states the bounded public adaptations and exclusions.

Author code is MIT, Copyright (c) 2026 Zeyu Fu (`LICENSE`). Authored manuscript
text, authored figures, derived result tables and documentation are CC BY 4.0
(`LICENSE-CONTENT.md`). These licenses do not relicense CPTAC/UMICH data, public
sample/protein identifiers, third-party packages, journal/TeX classes, fonts or
other upstream material. See `THIRD_PARTY_NOTICES.md` and `DATA_SOURCES.md`.

Private submission/cover-letter/admin files, all backups, old manuscript faces,
obsolete website assets, unrelated studies and third-party raw proteomics are
excluded. Scientific HOLD remains after every successful command.
