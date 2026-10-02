#!/usr/bin/env bash
# Public default verifies stored evidence and renders supplied PDFs; no PDF edits.
# --build explicitly regenerates figures/tables and compiles Article/SI locally.
set -euo pipefail
cd "$(dirname "$0")/../.."
export PYTHONDONTWRITEBYTECODE=1
python reproducibility/verify_method_contract.py
if [[ "${1:-}" == "--build" ]]; then
  export SOURCE_DATE_EPOCH=0 FORCE_SOURCE_DATE=1 TZ=UTC
  python manuscript/scripts/export_revision.py
  python manuscript/scripts/write_lock_methods_tex.py
  Rscript manuscript/scripts/R/make_b28_figures.R
  python manuscript/scripts/make_study_design.py
  python manuscript/scripts/compile_current.py
elif [[ $# -ne 0 ]]; then
  echo "Usage: rebuild_current.sh [--build]" >&2
  exit 2
fi
python analysis/revision_20260922/verify_candidate.py
python analysis/revision_20260922/verify_postprocess.py
python analysis/revision_20260922/verify_side_studies.py
python manuscript/scripts/check_layout.py --render
python manuscript/scripts/verify_integration.py --reuse-numerics
python reproducibility/verify_method_contract.py
