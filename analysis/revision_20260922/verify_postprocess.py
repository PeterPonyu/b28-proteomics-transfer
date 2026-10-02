#!/usr/bin/env python3
"""Read-only, fold-by-fold checks of corrected B28 calibration and coverage."""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score

from verify_candidate import ece

OUT = Path(__file__).resolve().parent / "candidate-v1"


def read(name):
    with (OUT / name).open(newline="") as f:
        return list(csv.DictReader(f))


def metrics(y, p):
    return dict(n=len(y), n_tumor=int(y.sum()), n_normal=int(len(y)-y.sum()),
                auroc=float(roc_auc_score(y,p)), ece_10=float(ece(y,p)),
                brier=float(brier_score_loss(y,p)), log_loss=float(log_loss(y,p)),
                accuracy=float(np.mean((p>=.5)==y)))


def main():
    original = [r for r in read("predictions.csv") if r["feature_set"] == "accession" and
                r["protocol"] == "leave_one_cohort_out"]
    lookup = {(r["cohort"],r["sample_id"]):r for r in original}
    recal = read("candidate-recalibration.csv")
    summary = json.loads((OUT / "candidate-recalibration-summary.json").read_text())
    if len(recal) != len(lookup) or len(set((r["cohort"],r["sample_id"]) for r in recal)) != len(lookup):
        raise AssertionError("Wrong calibration sample inventory")
    for r in recal:
        old = lookup[(r["cohort"],r["sample_id"])]
        if (r["patient_id"],r["y"],float(r["p_before"])) != (
                old["patient_id"],old["y"],float(old["p"])):
            raise AssertionError("Recalibration does not use the corrected frozen outer prediction")
    for fold in summary["folds"]:
        selected = [r for r in recal if r["cohort"] == fold["held_out"]]
        y = np.array([int(r["y"]) for r in selected])
        if fold["n_calibration"] != len(recal)-len(selected) or fold["n_training_cohorts_inner"] != [2,3]:
            raise AssertionError("Calibration inner study mixture differs")
        for name,col in (("before","p_before"),("temperature","p_temperature"),("platt","p_platt")):
            p = np.array([float(r[col]) for r in selected])
            for key,val in metrics(y,p).items():
                if abs(val-fold["metrics"][name][key])>1e-12:
                    raise AssertionError(f"{fold['held_out']} {name} {key} mismatch")

    rows = read("sample-identity-and-folds.csv")
    y = np.array([int(r["y"]) for r in rows])
    p = np.array([float(lookup[(r["cohort"],r["sample_id"])]["p"]) for r in rows])
    cohort = np.array([r["cohort"] for r in rows])
    for cov in read("candidate-riskcoverage.csv"):
        order = np.lexsort((np.arange(len(p)),-np.abs(p-.5)))
        idx = order[:max(2,int(float(cov["coverage"])*len(y)))]
        if int(cov["n_retained"]) != len(idx):
            raise AssertionError("Retained count disagrees")
        if abs(float(cov["accuracy_retained"])-np.mean((p[idx]>=.5)==y[idx]))>1e-12:
            raise AssertionError("Accuracy at coverage disagrees")
        if abs(float(cov["accuracy_full"])-np.mean((p>=.5)==y))>1e-12:
            raise AssertionError("Full accuracy disagrees")
        counts=Counter(cohort[idx])
        if any(int(cov[f"n_{c}"])!=counts[c] for c in set(cohort)):
            raise AssertionError("Retained cohort composition disagrees")
    print("PASS: corrected outer prediction identity, all calibration folds/metrics, and all coverage cells")


if __name__ == "__main__":
    main()
