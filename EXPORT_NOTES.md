# Public export adaptations and exclusions (2026-10-02)

This is an explicit dependency-aware allowlist export, not a capsule mirror.
`provenance/export-source-ledger.json` pins the source working-tree bytes for
each selected input; it does not claim that Git HEAD contains those bytes.

Bounded public adaptations, made only to the export:

- Lazy CPTAC import permits network-free function/identity tests; the cohort
  class names/order and scientific logic remain the same. Full fitting requests
  use local `no_internet=True` source loading and may select `B28_FIT_ROOT` for
  separate, non-overwriting outputs. CPTAC import initialization itself may
  still perform a network lookup, so full fitting is not called offline.
- Stored-prediction validators default to included evidence; optional explicit
  `--upstream-root` checks exact source hashes. Fresh QA goes under `output/`,
  not over frozen original QA/results. Upstream verification is explicitly false
  unless that separate option succeeds.
- The existing integration validator checks included numerical/compiled inputs
  and freshly rendered output. It uses the historical OV ledger only as labelled
  context; it no longer requires an excluded private submission PDF or original
  machine-specific compilation logs. Newly generated build logs are checked
  when present. No test was skipped to hide an omitted dependency.
- The default known rebuild/build entry points verify and render supplied PDFs
  without rewriting them. `--build` is an explicit offline figure/TeX rebuild.
  No automatic contract rebinding or release writer can bless changed inputs.
- Public integrity tests and manifests refer to actual included files. The
  default full lock fails on missing/modified/unexpected payload, historical
  schemas, path traversal, symlinks, and changed scientific evidence. The
  original manifests are historical, not falsely declared unchanged-build
  proofs for the trimmed export.

Original prediction CSVs, corrected accession/gene/sample/patient/fold maps,
recorded numerical outputs, current Article/SI source and PDFs, consumed figures
and table fragments are copied byte-for-byte. Changes above affect path/loading,
verification scope/output and public inventory only; no model predictions,
scientific score, fitting algorithm or claim was added or changed.

Excluded: all root and manuscript `submission/` trees; covers and administrative
checklists; `.git`, `.omc`, caches and private source pins; all backups, old
manuscript faces and candidate summary PDFs/figures; old website assets; historical
gene-indexed programs/results/splits and obsolete diagnostics; non-B28 content;
third-party abundance, mapping and annotation data files; standalone commercial
fonts; ZIPs and other submission/backup archives. The only retained old result
context is the clearly labelled OV summary and old contract/QA JSON necessary to
explain the provenance boundary.

Scientific HOLD is unchanged. Source/figure/hash verification is not scientific
clearance; OV is not prediction-replayed; cancer/study is confounded; the shared
panel uses all-six-cohort feature availability; fixed-prediction ranges are not
future-study confidence intervals. New reserved Zenodo DOI
10.5281/zenodo.23101188 is unpublished until the main release workflow publishes
it. This task performs no network mutation, commit, push, ZIP or PDF edit.
