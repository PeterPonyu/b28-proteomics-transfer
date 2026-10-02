# Public current reproducibility

Start with root `README.md`: stored prediction/metric checks, offline supplied
PDF/figure verification and full fitting have different dependencies and claims.
The public current method contract is `method-contract.json`. It pins the actual
included tree; `verify_method_contract.py` checks it by default and rejects the
historical capsule schema. Root `SHA256SUMS.txt` independently covers payload and
control files (except itself).

The current partition is `splits/patient-grouped-corrected-20260922.csv`;
`analysis/fold_identity.py` replays it by cohort/sample identity and fails on
duplicates, missing labels, changed patients or cross-fold patient assignments.
The invalid 723-key historical split is not included. True-patient counts and
HNSCC cored-normal identities are tested against actual stored rows.

All current scientific programs and stored outputs are under
`../analysis/revision_20260922/`. The numerical QA retains 301 checked metric
cells with 39 pinned included dependencies. Stored predictions and their original
QA are frozen; fresh validation outputs live under `../output/verification/`.

Original capsule manifest/QA context is explicitly historical under
`../provenance/historical/`. Its excluded dependency references do not validate
the public export. No original-matrix rerun, old diagnostic promotion, fresh fit,
OV prediction replay, future-study interval or scientific-readiness attestation
is implied by successful reproduction.
