#!/usr/bin/env python3
"""Candidate B28 LOGO recalibration and selective coverage on the corrected matrix.

The original two-group split-half inner protocol is retained here as a
comparison of the same *procedure* on corrected sample/feature identities.
It does not eliminate the inner-vs-outer training cohort-count mismatch.
"""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import xgboost as xgb

import reanalyse as core

HERE = core.OUT


def read(name):
    with (HERE / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def write(name, data):
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def logit(p):
    p = np.clip(p, 1e-7, 1-1e-7)
    return np.log(p/(1-p))


def sigmoid(x):
    return np.where(x >= 0, 1/(1+np.exp(-x)), np.exp(x)/(1+np.exp(x)))


def fit_inner(X, y, tr, te):
    med = np.nanmedian(X[tr], axis=0)
    med = np.where(np.isnan(med), 0., med)
    train = np.where(np.isnan(X[tr]), med, X[tr]).astype(np.float32)
    test = np.where(np.isnan(X[te]), med, X[te]).astype(np.float32)
    clf = xgb.XGBClassifier(**core.XGB)
    clf.fit(train, y[tr])
    return clf.predict_proba(test)[:, 1]


def fit_temperature(z, y):
    def nll(t):
        p = np.clip(sigmoid(z/t), 1e-7, 1-1e-7)
        return float(-np.mean(y*np.log(p)+(1-y)*np.log(1-p)))
    return float(minimize_scalar(nll, bounds=(.05, 20.), method="bounded").x)


def main():
    # Both outputs must be new. The frozen original pipeline is never touched.
    for name in ("candidate-recalibration.csv", "candidate-riskcoverage.csv"):
        if (HERE / name).exists():
            raise SystemExit(f"Refusing to overwrite {HERE/name}")
    _, rows, y, group, _, _, _, _, X, _ = core.load()
    key_to_i = {(r["cohort"], r["sample_id"]): i for i, r in enumerate(rows)}
    old = read("predictions.csv")
    preds = {(r["cohort"], r["sample_id"]): r for r in old if
             r["feature_set"] == "accession" and r["protocol"] == "leave_one_cohort_out"}
    if len(preds) != len(y) or set(preds) != set(key_to_i):
        raise AssertionError("Corrected LOGO predictions missing or repeated")
    source = json.loads((HERE / "candidate-results.json").read_text())["results"]["accession"]["leave_one_cohort_out"]
    ordered = np.array([float(preds[(r["cohort"],r["sample_id"])]["p"]) for r in rows])
    if abs(roc_auc_score(y, ordered)-source["pooled"]["auroc"]) > 1e-12:
        raise AssertionError("Not the same outer predictions")

    risk = []
    tie_order = np.lexsort((np.arange(len(y)), -np.abs(ordered-.5)))
    for coverage in np.linspace(.1, 1., 10):
        idx = tie_order[:max(2, int(coverage*len(y)))]
        correct = (ordered[idx] >= .5) == y[idx]
        counts = Counter(group[idx])
        record = dict(coverage=float(coverage), n_retained=len(idx),
                      accuracy_retained=float(correct.mean()),
                      accuracy_full=float(np.mean((ordered >= .5) == y)))
        record.update({f"n_{cohort}":counts[cohort] for cohort in core.COHORTS})
        risk.append(record)
    write("candidate-riskcoverage.csv", risk)

    recal = []
    fold_details = []
    for held in core.COHORTS:
        train_groups = [c for c in core.COHORTS if c != held]
        group_a, group_b = train_groups[:2], train_groups[2:]
        cal_p, cal_y = [], []
        for fit_groups, eval_groups in ((group_a,group_b),(group_b,group_a)):
            tr = np.flatnonzero(np.isin(group, fit_groups))
            te = np.flatnonzero(np.isin(group, eval_groups))
            cal_p.append(fit_inner(X, y, tr, te))
            cal_y.append(y[te])
        calibration_score = logit(np.concatenate(cal_p))
        calibration_label = np.concatenate(cal_y)
        temperature = fit_temperature(calibration_score, calibration_label)
        platt = LogisticRegression(C=1e6, max_iter=2000)
        platt.fit(calibration_score.reshape(-1,1), calibration_label)
        te = np.flatnonzero(group == held)
        z = logit(ordered[te])
        p_t = sigmoid(z/temperature)
        p_platt = platt.predict_proba(z.reshape(-1,1))[:,1]
        for j, t, pl in zip(te,p_t,p_platt):
            r = rows[j]
            recal.append(dict(cohort=r["cohort"], sample_id=r["sample_id"],
                              patient_id=r["patient_id"], y=int(y[j]),
                              p_before=float(ordered[j]), p_temperature=float(t),
                              p_platt=float(pl), temperature=temperature))
        fold_details.append(dict(held_out=held, n_calibration=len(calibration_label),
                                 n_training_cohorts_inner=[2,3], n_training_cohorts_outer=5,
                                 fitted_temperature=temperature))
        print(f"recalibrated held-out {held}, n={len(te)}, T={temperature:.3f}", flush=True)
    write("candidate-recalibration.csv", recal)
    for row in fold_details:
        cohort = row["held_out"]
        subset = [r for r in recal if r["cohort"] == cohort]
        lab = np.array([r["y"] for r in subset])
        row["metrics"] = {name: core.metrics(lab, np.array([r[field] for r in subset]))
                          for name, field in (("before","p_before"),
                                              ("temperature","p_temperature"),("platt","p_platt"))}
    (HERE / "candidate-recalibration-summary.json").write_text(json.dumps(dict(
        status="candidate, original split-half calibration procedure rerun on corrected objects; manuscript HOLD",
        caveat="inner 2/3-cohort fit differs from 5-cohort outer model, so not a universal recalibration mechanism test",
        folds=fold_details), indent=2) + "\n")
    print("saved candidate recalibration and coverage", flush=True)


if __name__ == "__main__":
    main()
