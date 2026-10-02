#!/usr/bin/env python3
"""Paired true-patient resampling of the three B28 *stored* prediction sets.

These percentile ranges condition on fitted models and six selected cohorts;
they are not confidence intervals for refitting the full algorithm or for
new independent studies. A patient is resampled within their source cohort.
"""
from __future__ import annotations

import csv
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

HERE = Path(os.environ.get('B28_FIT_ROOT', str(Path(__file__).resolve().parent))) / "candidate-v1"
PROTOCOLS = ("sample_random", "patient_grouped", "leave_one_cohort_out")
SEED = 20260922
N_DRAW = 2000


def read(name):
    with (HERE / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def draw_metrics(indices, ys, scores, cohorts):
    by_cohort = []
    for c in sorted(set(cohorts)):
        subset = indices[cohorts[indices] == c]
        if len(set(ys[subset])) < 2:
            return None
        by_cohort.append([roc_auc_score(ys[subset], scores[p][subset]) for p in PROTOCOLS])
    per = np.array(by_cohort)
    macro = per.mean(axis=0)
    pooled = np.array([roc_auc_score(ys[indices], scores[p][indices]) for p in PROTOCOLS])
    return macro, pooled


def main():
    target = HERE / "fixed-prediction-patient-resampling.json"
    if target.exists():
        raise SystemExit(f"Refusing to overwrite: {target}")
    rows = read("sample-identity-and-folds.csv")
    records = read("predictions.csv")
    keys = [(r["cohort"], r["sample_id"]) for r in rows]
    scores = {}
    for protocol in PROTOCOLS:
        mapping = {(r["cohort"], r["sample_id"]): float(r["p"]) for r in records
                   if r["feature_set"] == "accession" and r["protocol"] == protocol}
        if len(mapping) != len(rows) or set(mapping) != set(keys):
            raise AssertionError(f"Not paired on full inventory: {protocol}")
        scores[protocol] = np.array([mapping[k] for k in keys])
    ys = np.array([int(r["y"]) for r in rows])
    cohorts = np.array([r["cohort"] for r in rows])
    patients = np.array([r["patient_id"] for r in rows])
    by_cohort_patient = defaultdict(lambda: defaultdict(list))
    for i, (c, patient) in enumerate(zip(cohorts, patients)):
        by_cohort_patient[c][patient].append(i)
    patient_lists = {c: list(d.values()) for c, d in by_cohort_patient.items()}
    point = draw_metrics(np.arange(len(rows)), ys, scores, cohorts)
    if point is None:
        raise AssertionError("Original cohort has one class")
    rng = np.random.default_rng(SEED)
    draws = []
    for _ in range(N_DRAW):
        idx = np.concatenate([np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
                              for c, groups in patient_lists.items()])
        metric = draw_metrics(idx, ys, scores, cohorts)
        if metric is not None:
            draws.append(metric)
    if len(draws) != N_DRAW:
        raise AssertionError(f"Unexpected single-class bootstrap draw: {len(draws)}/{N_DRAW}")
    result = {}
    for name, j in (("macro_equal_cohort",0),("pooled",1)):
        original = point[j]
        samples = np.array([v[j] for v in draws])
        contrasts = {
            "sample_minus_logo": (0,2), "patient_minus_logo": (1,2),
            "sample_minus_patient": (0,1),
        }
        result[name] = dict(
            auroc={p:dict(point=float(original[i]), conditional_percentile_range=[float(v) for v in np.percentile(samples[:,i], [2.5,97.5])])
                   for i, p in enumerate(PROTOCOLS)},
            contrasts={label:dict(point=float(original[a]-original[b]),
                                  conditional_percentile_range=[float(v) for v in np.percentile(samples[:,a]-samples[:,b], [2.5,97.5])])
                       for label,(a,b) in contrasts.items()})
    target.write_text(json.dumps(dict(status="fixed-prediction conditional sensitivity, not a model-refit or new-study CI",
                                      n_draws=N_DRAW, seed=SEED, unit="true patient nested within fixed source cohort",
                                      results=result), indent=2)+"\n")
    print("saved", target)


if __name__ == "__main__":
    main()
