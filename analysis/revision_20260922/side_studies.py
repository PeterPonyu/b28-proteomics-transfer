#!/usr/bin/env python3
"""B28 candidate LOGO sensitivities on the corrected stable-ENSP matrix.

This is a separate, non-overwriting candidate run. The full-matrix median arm
deliberately leaks held-out feature values; the L2 arm uses train-only median
imputation and train-only standardization. Neither arm changes the source
identity, feature inventory, or held-out fold of candidate-v1.
"""
from __future__ import annotations

import csv
import hashlib
import json
import warnings
from importlib.metadata import version
from pathlib import Path

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

import reanalyse as core

HERE = Path(__file__).resolve().parent
BASE = core.OUT
OUT = core.FIT_ROOT / "candidate-v2-side-studies"
ARMS = ("full_data_median_leak_xgb", "train_median_l2")


def fill_with_median(X: np.ndarray, train: np.ndarray, test: np.ndarray,
                     *, leak: bool) -> tuple[np.ndarray, np.ndarray]:
    reference = X if leak else X[train]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-missing columns -> zero below
        med = np.nanmedian(reference, axis=0)
    med = np.where(np.isnan(med), 0., med)
    return (np.where(np.isnan(X[train]), med, X[train]).astype(np.float32),
            np.where(np.isnan(X[test]), med, X[test]).astype(np.float32))


def predict_arm(arm: str, X: np.ndarray, y: np.ndarray,
                train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, dict]:
    if arm not in ARMS:
        raise ValueError(f"Unknown arm: {arm}")
    Xtr, Xte = fill_with_median(X, train, test,
                                leak=arm == "full_data_median_leak_xgb")
    if arm == "full_data_median_leak_xgb":
        clf = xgb.XGBClassifier(**core.XGB)
        clf.fit(Xtr, y[train])
        return clf.predict_proba(Xte)[:, 1], {"converged": None}
    scale = StandardScaler(with_mean=True, with_std=True)
    Xtr = scale.fit_transform(Xtr)
    Xte = scale.transform(Xte)
    clf = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=500)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        clf.fit(Xtr, y[train])
    converged = not any(issubclass(w.category, ConvergenceWarning) for w in caught)
    return clf.predict_proba(Xte)[:, 1], {
        "converged": converged,
        "iterations": int(clf.n_iter_[0]),
        "max_abs_train_center": float(np.max(np.abs(Xtr.mean(axis=0)))),
    }


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite existing candidate: {OUT}")
    baseline = json.loads((BASE / "candidate-results.json").read_text())
    manifest = read_csv(BASE / "sample-identity-and-folds.csv")
    frozen = read_csv(BASE / "predictions.csv")
    _, rows, y, groups, patients, accessions, _, _, X, _ = core.load()
    if len(accessions) != 8906 or len(rows) != 1191 or len(rows) != len(manifest):
        raise AssertionError("Corrected cohort object changed")
    keys = [(r["cohort"], r["sample_id"]) for r in rows]
    if keys != [(r["cohort"], r["sample_id"]) for r in manifest]:
        raise AssertionError("Source order and candidate-v1 sample identity differ")
    if any((int(m["y"]), m["patient_id"]) != (int(y[i]), patients[i])
           for i, m in enumerate(manifest)):
        raise AssertionError("Label or patient key changed")
    source = {(r["cohort"], r["sample_id"]): r for r in frozen if
              r["feature_set"] == "accession" and r["protocol"] == "leave_one_cohort_out"}
    if len(source) != len(rows) or set(source) != set(keys):
        raise AssertionError("Missing/duplicate reference LOGO predictions")
    splits, assignments = core.split_manifest(y, groups, patients, rows)
    if any(int(source[k]["fold"]) != assignments[("leave_one_cohort_out", i)]
           for i, k in enumerate(keys)):
        raise AssertionError("Reference folds differ from corrected cohort object")
    if any(source[k]["patient_id"] != patients[i] or int(source[k]["y"]) != y[i]
           for i, k in enumerate(keys)):
        raise AssertionError("Reference identities/labels differ")

    records, fit_details = [], []
    for arm in ARMS:
        for fold, (tr, te) in enumerate(splits["leave_one_cohort_out"], 1):
            prob, fit = predict_arm(arm, X, y, tr, te)
            held = str(groups[te[0]])
            if len(set(groups[te])) != 1 or not np.isfinite(prob).all():
                raise AssertionError("Invalid held-out cohort or predictions")
            if arm == "train_median_l2" and not fit["converged"]:
                raise RuntimeError(f"L2 did not converge for {held}: {fit}")
            fit_details.append(dict(arm=arm, fold=fold, held_out=held,
                                    n_train=len(tr), n_test=len(te), **fit))
            for i, p in zip(te, prob):
                records.append(dict(arm=arm, protocol="leave_one_cohort_out",
                                    fold=fold, cohort=rows[i]["cohort"],
                                    sample_id=rows[i]["sample_id"],
                                    patient_id=rows[i]["patient_id"],
                                    tissue=rows[i]["tissue"], y=int(y[i]), p=float(p)))
            print(f"{arm}: {held}, n={len(te)}, AUROC={core.metrics(y[te], prob)['auroc']:.4f}", flush=True)

    results = core.summarize([dict(feature_set=r["arm"], **{k: v for k, v in r.items() if k != "arm"})
                              for r in records], rows)
    results = {arm: results[arm]["leave_one_cohort_out"] for arm in ARMS}
    ref = baseline["results"]["accession"]["leave_one_cohort_out"]
    contrasts = {}
    for arm in ARMS:
        contrasts[arm] = {place: {metric: results[arm][place][metric]-ref[place][metric]
                                   for metric in ("auroc", "ece_10", "brier", "log_loss")}
                          for place in ("macro_equal_cohort", "pooled")}
    OUT.mkdir(parents=True)
    core.write_csv(OUT / "predictions.csv", list(records[0]), records)
    report = dict(status="candidate side studies only; original manuscript HOLD",
                  input_sha256=baseline["input_sha256"],
                  candidate_v1_sha256={name: hashlib.sha256((BASE/name).read_bytes()).hexdigest()
                                       for name in ("candidate-results.json", "predictions.csv",
                                                    "sample-identity-and-folds.csv")},
                  shared_accessions=len(accessions), n_samples=len(rows),
                  protocol="same six corrected LOGO folds as candidate-v1",
                  full_data_median_leak="deliberately includes all held-out feature values, no held-out labels",
                  train_median_l2=dict(penalty="l2", C=1.0, solver="lbfgs", max_iter=500,
                                       standardization="fit on each training fold only"),
                  software={name: version(name) for name in
                            ("cptac", "xgboost", "numpy", "pandas", "scikit-learn")},
                  fit_details=fit_details, results=results, contrast_against_accession_xgb=contrasts)
    (OUT / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print("saved", OUT, flush=True)


if __name__ == "__main__":
    main()
