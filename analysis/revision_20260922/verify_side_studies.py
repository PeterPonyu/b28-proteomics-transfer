#!/usr/bin/env python3
"""Independently recompute B28 candidate-v2 metrics from row-level outputs."""
from __future__ import annotations

import csv
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score

HERE = Path(__file__).resolve().parent
BASE = HERE / "candidate-v1"
OUT = HERE / "candidate-v2-side-studies"
COHORTS = ("Ccrcc", "Luad", "Ucec", "Hnscc", "Lscc", "Pdac")


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def evaluate(rows):
    y = np.array([int(r["y"]) for r in rows])
    p = np.array([float(r["p"]) for r in rows])
    edges = np.linspace(0, 1, 11)
    ece = 0.
    for k in range(10):
        membership = (p >= edges[k]) & ((p < edges[k+1]) if k < 9 else (p <= edges[k+1]))
        if np.any(membership):
            ece += float(np.mean(membership) * abs(p[membership].mean()-y[membership].mean()))
    return dict(auroc=float(roc_auc_score(y, p)), ece_10=ece,
                brier=float(brier_score_loss(y, p)), log_loss=float(log_loss(y, p)))


def check(saved, rows):
    actual = evaluate(rows)
    for key, value in actual.items():
        if abs(saved[key]-value) > 1e-12:
            raise AssertionError(f"{key}: {saved[key]} vs independent {value}")
    return actual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream-root', type=Path,
                        help='Explicit cptac/data root for optional upstream hash verification')
    args = parser.parse_args()
    report = json.loads((OUT / "results.json").read_text())
    ref = json.loads((BASE / "candidate-results.json").read_text())
    identity = read(BASE / "sample-identity-and-folds.csv")
    source = read(BASE / "predictions.csv")
    pred = read(OUT / "predictions.csv")
    lookup = {(r["cohort"], r["sample_id"]): r for r in identity}
    base_lookup = {(r["cohort"], r["sample_id"]): r for r in source if
                   r["feature_set"] == "accession" and r["protocol"] == "leave_one_cohort_out"}
    if set(lookup) != set(base_lookup) or len(lookup) != 1191:
        raise AssertionError("Incorrect reference population")
    arms = defaultdict(list)
    for r in pred:
        key = (r["cohort"], r["sample_id"])
        if key not in lookup:
            raise AssertionError("Unknown sample")
        orig, prior = lookup[key], base_lookup[key]
        if (r["protocol"], r["patient_id"], r["tissue"], r["y"], r["fold"]) != (
                "leave_one_cohort_out", orig["patient_id"], orig["tissue"], orig["y"],
                orig["fold_leave_one_cohort_out"]):
            raise AssertionError(f"Wrong prediction identity/fold for {key}")
        if (r["y"], r["fold"]) != (prior["y"], prior["fold"]):
            raise AssertionError("Not comparable with corrected reference")
        if not 0 <= float(r["p"]) <= 1:
            raise AssertionError("Invalid probability")
        arms[r["arm"]].append(r)
    if set(arms) != {"full_data_median_leak_xgb", "train_median_l2"}:
        raise AssertionError("Missing/unknown arm")
    for arm, records in arms.items():
        if Counter((r["cohort"], r["sample_id"]) for r in records) != Counter(lookup.keys()):
            raise AssertionError("Arm did not predict every sample exactly once")
        result = report["results"][arm]
        check(result["pooled"], records)
        checks = []
        for cohort in COHORTS:
            subset = [r for r in records if r["cohort"] == cohort]
            checks.append(check(result["by_cohort"][cohort], subset))
        for metric in ("auroc", "ece_10", "brier", "log_loss"):
            macro = float(np.mean([m[metric] for m in checks]))
            if abs(result["macro_equal_cohort"][metric]-macro) > 1e-12:
                raise AssertionError("Invalid equal-cohort macro")
            for place in ("macro_equal_cohort", "pooled"):
                target = result[place][metric]-ref["results"]["accession"]["leave_one_cohort_out"][place][metric]
                saved = report["contrast_against_accession_xgb"][arm][place][metric]
                if abs(target-saved) > 1e-12:
                    raise AssertionError("Incorrect contrast")
    if len(report["fit_details"]) != 12 or any(
            item["converged"] is False for item in report["fit_details"]):
        raise AssertionError("Missing fold or nonconverged linear model")
    if args.upstream_root is not None:
        for path, sha in report["input_sha256"].items():
            local = args.upstream_root / Path(path).parent.name / Path(path).name
            if not local.is_file() or hashlib.sha256(local.read_bytes()).hexdigest() != sha:
                raise AssertionError(f"Upstream source absent or changed: {local}")
    for name, sha in report["candidate_v1_sha256"].items():
        if hashlib.sha256((BASE/name).read_bytes()).hexdigest() != sha:
            raise AssertionError("Candidate reference changed")
    qa = dict(status="candidate-v2 row-level verification passed; manuscript HOLD",
              n_predictions=len(pred), n_arms=len(arms), n_folds=len(report["fit_details"]),
              upstream_verified_in_this_run=args.upstream_root is not None,
              output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in (OUT / "predictions.csv", OUT / "results.json")})
    output = HERE.parents[1] / 'output/verification'
    output.mkdir(parents=True, exist_ok=True)
    (output / "side-studies-QA.json").write_text(json.dumps(qa, indent=2) + "\n")
    print(json.dumps(qa, indent=2))


if __name__ == "__main__":
    main()
