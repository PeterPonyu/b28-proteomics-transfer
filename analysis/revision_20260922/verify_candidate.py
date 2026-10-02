#!/usr/bin/env python3
"""Independent checks on saved B28 candidate files, without fitting models."""
from __future__ import annotations

import csv
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score

HERE = Path(__file__).resolve().parent / "candidate-v1"
COHORTS = {"Ccrcc", "Luad", "Ucec", "Hnscc", "Lscc", "Pdac"}


def read(name):
    with (HERE / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def ece(y, p):
    edges = np.linspace(0, 1, 11)
    total = 0.
    for i in range(10):
        mask = (p >= edges[i]) & (p < edges[i+1] if i < 9 else p <= edges[i+1])
        if mask.any():
            total += mask.sum()/len(p)*abs(np.mean(y[mask])-np.mean(p[mask]))
    return total


def compare(saved, y, p):
    values = dict(auroc=roc_auc_score(y, p), ece_10=ece(y, p),
                  brier=brier_score_loss(y, p), log_loss=log_loss(y, p))
    for metric, value in values.items():
        if abs(saved[metric]-value) > 1e-12:
            raise AssertionError(f"{metric}: stored {saved[metric]}, recomputed {value}")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream-root', type=Path,
                        help='Explicit cptac/data root for optional upstream hash verification')
    args = parser.parse_args()
    result = json.loads((HERE / "candidate-results.json").read_text())
    rows = read("sample-identity-and-folds.csv")
    pred = read("predictions.csv")
    accession = read("primary-accession-map.csv")
    gene = read("gene-accession-audit.csv")
    keys = [(r["cohort"], r["sample_id"]) for r in rows]
    if len(rows) != 1191 or len(set(keys)) != len(rows):
        raise AssertionError("Sample identity inventory incomplete or nonunique")
    keyed = dict(zip(keys, rows))
    for r in rows:
        suffix = r["sample_id"].rsplit(".", 1)
        if r["cohort"] not in COHORTS:
            raise AssertionError("Unknown cohort")
        expected = 0 if len(suffix) == 2 and suffix[1] in {"N", "C"} else 1
        tissue = "tumor" if expected else ("cored_normal" if suffix[1] == "C" else "adjacent_normal")
        if suffix[-1] == "C" and r["cohort"] != "Hnscc":
            raise AssertionError("Cored normal outside HNSCC")
        if (int(r["y"]), r["tissue"], r["patient_id"]) != (
            expected, tissue, f'{r["cohort"]}|{suffix[0] if not expected else r["sample_id"]}'):
            raise AssertionError(f"Incorrect sample identity: {r}")
    if Counter(r["tissue"] for r in rows) != Counter(tumor=689, adjacent_normal=496, cored_normal=6):
        raise AssertionError("Tissue counts changed")
    by_patient = defaultdict(list)
    for r in rows:
        by_patient[r["patient_id"]].append(r)
    if len(by_patient) != 719:
        raise AssertionError("Patient count changed")
    for protocol in ("patient_grouped", "leave_one_cohort_out"):
        if any(len({r[f"fold_{protocol}"] for r in rr}) != 1 for rr in by_patient.values()):
            raise AssertionError(f"A true patient was split: {protocol}")
    if Counter(r["tissue"] for r in rows if r["cohort"] == "Hnscc") != Counter(
            tumor=110, adjacent_normal=66, cored_normal=6):
        raise AssertionError("HNSCC annotation changed")

    acc_by = defaultdict(dict)
    for r in accession:
        acc_by[r["stable_accession"]][r["cohort"]] = r
        if r["original_accession"].split(".")[0] != r["stable_accession"]:
            raise AssertionError("Stable/full accession mismatch")
    if len(acc_by) != 8906 or any(set(rr) != COHORTS or
                                  len({v["gene"] for v in rr.values()}) != 1
                                  for rr in acc_by.values()):
        raise AssertionError("Not a shared protein identity")
    first = defaultdict(dict)
    for r in gene:
        first[r["gene"]][r["cohort"]] = r
    if len(first) != 9675 or any(set(rr) != COHORTS for rr in first.values()):
        raise AssertionError("Old gene inventory differs")
    ambiguous = sum(len({r["first_stable_accession"] for r in rr.values()}) > 1
                    for rr in first.values())
    if ambiguous != 1234:
        raise AssertionError(f"Old mismatched-gene count is {ambiguous}")

    prediction_groups = defaultdict(list)
    for r in pred:
        key = (r["cohort"], r["sample_id"])
        if key not in keyed:
            raise AssertionError("Unknown prediction row")
        original = keyed[key]
        if (r["patient_id"], r["tissue"], r["y"], r["fold"]) != (
                original["patient_id"], original["tissue"], original["y"],
                original[f'fold_{r["protocol"]}']):
            raise AssertionError("Prediction metadata/fold differs from manifest")
        if not (0 <= float(r["p"]) <= 1):
            raise AssertionError("Probability outside [0,1]")
        prediction_groups[(r["feature_set"], r["protocol"])].append(r)
    for (feature, protocol), rr in prediction_groups.items():
        if len(rr) != len(rows) or Counter((r["cohort"], r["sample_id"]) for r in rr) != Counter(keys):
            raise AssertionError(f"Output is not one-to-one: {feature} {protocol}")
        summary = result["results"][feature][protocol]
        y = np.array([int(r["y"]) for r in rr]); p = np.array([float(r["p"]) for r in rr])
        compare(summary["pooled"], y, p)
        for cohort in COHORTS:
            mask = np.array([r["cohort"] == cohort for r in rr])
            compare(summary["by_cohort"][cohort], y[mask], p[mask])
        macro = np.mean([summary["by_cohort"][c]["auroc"] for c in COHORTS])
        if abs(macro-summary["macro_equal_cohort"]["auroc"]) > 1e-12:
            raise AssertionError("Macro AUROC not equal-cohort")
    if set(prediction_groups) != ({("accession", p) for p in (
            "sample_random", "patient_grouped", "leave_one_cohort_out")} |
            {(s, "leave_one_cohort_out") for s in (
            "consistent_gene_first_sensitivity", "global_summaries", "missingness")}):
        raise AssertionError("A planned prediction set is absent")

    source_digests = result["input_sha256"]
    if args.upstream_root is not None:
        for path, sha in source_digests.items():
            local = args.upstream_root / Path(path).parent.name / Path(path).name
            if not local.is_file() or hashlib.sha256(local.read_bytes()).hexdigest() != sha:
                raise AssertionError(f"Upstream source absent or changed: {local}")
    qa = dict(status="candidate identity, predictions and metrics independently verified; manuscript HOLD",
              n_predictions=len(pred), n_prediction_sets=len(prediction_groups),
              n_patients=len(by_patient), n_true_patient_split_violations=0,
              n_old_gene_first_accession_mismatches=ambiguous,
              input_sha256=source_digests,
              upstream_verified_in_this_run=args.upstream_root is not None,
              output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in HERE.glob("*.csv")})
    output = HERE.parents[2] / 'output/verification'
    output.mkdir(parents=True, exist_ok=True)
    (output / "candidate-QA.json").write_text(json.dumps(qa, indent=2) + "\n")
    print(json.dumps({k: v for k, v in qa.items() if k not in ("input_sha256", "output_sha256")}, indent=2))


if __name__ == "__main__":
    main()
