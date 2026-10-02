# Historical records only: not the public current dependency closure

These narrowly allowlisted JSON records are byte-identical historical context.
No private submission document, old manuscript PDF, cover letter, administrative
checklist or backup is included.

- `capsule-method-contract.json` and `capsule-release-current.json` are the
  original 2026-09-23 capsule contracts. They refer to intentionally excluded
  files (including backups and a historical submission PDF) and cannot verify
  this public tree. Use root `PUBLIC-MANIFEST.json` and the public
  `reproducibility/method-contract.json` instead.
- `capsule-integration-QA.json` is prior completed integration context. Its old
  file/raster assertions are not a fresh public-run attestation.
- `ov-summary-20260911.json` is the historical OV plex result referred to in SI
  S7. Its old headline references and old bootstrap endpoints are historical,
  not corrected six-cohort results or new inference. The OV prediction rows
  were not replayed. `capsule-ov-source-check.json` records the prior source/hash
  and annotation inventory check only. This export verifies ledger identity,
  not the omitted upstream raw files.

Absolute input paths in these records, and in the frozen result manifests, are
historical provenance strings. No default public verification dereferences them.
The exact hashes remain useful to match separately obtained upstream files.
