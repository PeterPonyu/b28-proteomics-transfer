#!/usr/bin/env python3
"""Candidate B28 reanalysis: true tissue/patient identity and common ENSP.

This is intentionally separate from the frozen manuscript pipeline. It reads
the six cached UMICH tables through cptac and refuses ambiguous identities.
Run with the recorded CPTAC environment; no legacy result is overwritten.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import warnings
from collections import defaultdict
from importlib.metadata import version
from pathlib import Path

import numpy as np
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut, StratifiedKFold
import xgboost as xgb

HERE = Path(__file__).resolve().parent
FIT_ROOT = Path(os.environ.get('B28_FIT_ROOT', str(HERE)))
OUT = FIT_ROOT / "candidate-v1"
COHORTS = {
    "Ccrcc": "Ccrcc", "Luad": "Luad", "Ucec": "Ucec",
    "Hnscc": "Hnscc", "Lscc": "Lscc", "Pdac": "Pdac",
}
XGB = dict(n_estimators=150, max_depth=4, learning_rate=.1, subsample=.8,
           colsample_bytree=.3, eval_metric="logloss", tree_method="hist",
           n_jobs=4, verbosity=0)


def stable_accession(raw: str) -> str:
    match = re.fullmatch(r"(ENSP\d+)(?:\.\d+)?", str(raw))
    if match is None:
        raise ValueError(f"Unsupported protein accession: {raw!r}")
    return match.group(1)


def tissue_patient(cohort: str, sample: str) -> tuple[str, int, str]:
    sid = str(sample)
    if sid.endswith(".N"):
        return "adjacent_normal", 0, f"{cohort}|{sid[:-2]}"
    if sid.endswith(".C") and cohort == "Hnscc":
        return "cored_normal", 0, f"{cohort}|{sid[:-2]}"
    if re.fullmatch(r"(?:C3[NL]-\d+|11LU\d+|NX\d+)", sid) is None:
        raise ValueError(f"Unknown sample status: {cohort} {sid!r}")
    return "tumor", 1, f"{cohort}|{sid}"


def ece(y: np.ndarray, p: np.ndarray) -> float:
    edges = np.linspace(0, 1, 11)
    total = 0.
    for i in range(10):
        mask = (p >= edges[i]) & (p < edges[i+1] if i < 9 else p <= edges[i+1])
        if mask.any():
            total += mask.mean() * abs(p[mask].mean() - y[mask].mean())
    return float(total)


def metrics(y: np.ndarray, p: np.ndarray) -> dict:
    if len(np.unique(y)) < 2:
        raise ValueError("A cohort/fold has only one class")
    return dict(n=int(len(y)), n_tumor=int(y.sum()), n_normal=int(len(y)-y.sum()),
                auroc=float(roc_auc_score(y, p)), ece_10=ece(y, p),
                brier=float(brier_score_loss(y, p)), log_loss=float(log_loss(y, p)),
                accuracy=float(np.mean((p >= .5) == y)))


def load():
    # cptac refreshes its remote index on import. Offline tests must not import
    # it; an explicitly requested full fit still requires the upstream package.
    import cptac

    frames = {}
    for cohort, source in COHORTS.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=Warning)
            frame = getattr(cptac, source)(no_internet=True).get_proteomics("umich")
        if frame.columns.nlevels != 2:
            raise ValueError(f"Expected gene/accession MultiIndex: {cohort}")
        if not frame.index.is_unique:
            raise ValueError(f"Duplicate sample: {cohort}")
        by_accession = {}
        by_gene = defaultdict(list)
        for i, (gene, raw) in enumerate(frame.columns):
            acc = stable_accession(raw)
            if acc in by_accession:
                raise ValueError(f"Duplicate stable accession: {cohort} {acc}")
            by_accession[acc] = (i, str(gene), str(raw))
            by_gene[str(gene)].append(acc)
        metadata = [dict(cohort=cohort, sample_id=str(sid),
                         tissue=tissue_patient(cohort, str(sid))[0],
                         y=tissue_patient(cohort, str(sid))[1],
                         patient_id=tissue_patient(cohort, str(sid))[2])
                    for sid in frame.index]
        frames[cohort] = dict(table=frame, by_accession=by_accession,
                              by_gene=by_gene, metadata=metadata)
        print(f"loaded {cohort}: {frame.shape}", flush=True)

    accessions = sorted(set.intersection(*(set(d["by_accession"]) for d in frames.values())))
    shared_genes = sorted(set.intersection(*(set(d["by_gene"]) for d in frames.values())))
    first_consistent = sorted(g for g in shared_genes if len({
        frames[c]["by_gene"][g][0] for c in COHORTS}) == 1)
    for acc in accessions:
        if len({frames[c]["by_accession"][acc][1] for c in COHORTS}) != 1:
            raise ValueError(f"One accession mapped to different gene names: {acc}")
    if (len(accessions), len(shared_genes), len(first_consistent)) != (8906, 9675, 8441):
        raise ValueError("Feature inventory differs from reviewed source; inspect before fitting")

    rows = [r for c in COHORTS for r in frames[c]["metadata"]]
    y = np.array([r["y"] for r in rows], dtype=int)
    groups = np.array([r["cohort"] for r in rows])
    patients = np.array([r["patient_id"] for r in rows])
    if (len(rows), int(y.sum()), len(set(patients))) != (1191, 689, 719):
        raise ValueError("Sample/tissue/patient inventory differs from review")
    X = np.vstack([d["table"].iloc[:, [d["by_accession"][a][0] for a in accessions]].to_numpy(dtype=float)
                   for d in frames.values()])
    X_consistent = np.vstack([d["table"].iloc[:, [d["by_accession"][d["by_gene"][g][0]][0]
                                                   for g in first_consistent]].to_numpy(dtype=float)
                               for d in frames.values()])
    if X.shape != (1191, 8906) or X_consistent.shape != (1191, 8441):
        raise AssertionError("Wrong feature matrix shape")
    return frames, rows, y, groups, patients, accessions, shared_genes, first_consistent, X, X_consistent


def split_manifest(y, groups, patients, rows):
    sample = list(StratifiedKFold(5, shuffle=True, random_state=28).split(np.zeros(len(y)), y))
    patient = list(GroupKFold(5).split(np.zeros(len(y)), y, patients))
    logo = list(LeaveOneGroupOut().split(np.zeros(len(y)), y, groups))
    splits = dict(sample_random=sample, patient_grouped=patient, leave_one_cohort_out=logo)
    assignments = {}
    for protocol, folds in splits.items():
        seen = set()
        for i, (tr, te) in enumerate(folds, 1):
            if set(tr) & set(te) or not set(tr) or not set(te):
                raise AssertionError("Invalid train/test fold")
            if protocol != "sample_random" and set(patients[tr]) & set(patients[te]):
                raise AssertionError("True patient crosses fold")
            for j in te:
                if j in seen:
                    raise AssertionError("Sample evaluated twice")
                seen.add(j)
                assignments[(protocol, int(j))] = i
        if seen != set(range(len(y))):
            raise AssertionError("Incomplete split")
    return splits, assignments


def predict(X, y, rows, splits, feature, records):
    for protocol, folds in splits.items():
        for fold, (tr, te) in enumerate(folds, 1):
            if feature != "accession" and protocol != "leave_one_cohort_out":
                continue
            if feature == "global_summaries":
                # Source values are UMICH reference-subtracted log-ratio measures.
                F = np.column_stack((np.nanmean(X, axis=1), np.nanmedian(X, axis=1),
                                     np.mean(~np.isnan(X), axis=1)))
            elif feature == "missingness":
                F = np.isnan(X).astype(np.float32)
            else:
                F = X
            med = np.nanmedian(F[tr], axis=0)
            med = np.where(np.isnan(med), 0., med)
            train = np.where(np.isnan(F[tr]), med, F[tr]).astype(np.float32)
            test = np.where(np.isnan(F[te]), med, F[te]).astype(np.float32)
            clf = xgb.XGBClassifier(**XGB)
            clf.fit(train, y[tr])
            prob = clf.predict_proba(test)[:, 1]
            for j, p in zip(te, prob):
                records.append(dict(feature_set=feature, protocol=protocol, fold=fold,
                                    cohort=rows[j]["cohort"], sample_id=rows[j]["sample_id"],
                                    patient_id=rows[j]["patient_id"], tissue=rows[j]["tissue"],
                                    y=int(y[j]), p=float(p)))
            print(f"{feature} {protocol} {fold}: {len(te)} held out", flush=True)


def summarize(records, rows):
    summary = {}
    source_ids = {(r["cohort"], r["sample_id"]) for r in rows}
    for feature in sorted({r["feature_set"] for r in records}):
        summary[feature] = {}
        for protocol in sorted({r["protocol"] for r in records if r["feature_set"] == feature}):
            item = [r for r in records if r["feature_set"] == feature and r["protocol"] == protocol]
            if len(item) != len(rows) or {(r["cohort"], r["sample_id"]) for r in item} != source_ids:
                raise AssertionError("Output does not cover each sample exactly once")
            by_cohort = {c: metrics(np.array([r["y"] for r in item if r["cohort"] == c]),
                                    np.array([r["p"] for r in item if r["cohort"] == c]))
                         for c in COHORTS}
            pooled = metrics(np.array([r["y"] for r in item]), np.array([r["p"] for r in item]))
            macro = {key: float(np.mean([m[key] for m in by_cohort.values()]))
                     for key in ("auroc", "ece_10", "brier", "log_loss", "accuracy")}
            summary[feature][protocol] = dict(by_cohort=by_cohort, macro_equal_cohort=macro,
                                               pooled=pooled)
    return summary


def write_csv(path, fields, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    import cptac

    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite existing candidate output: {OUT}")
    OUT.mkdir(parents=True)
    frames, rows, y, groups, patients, accessions, genes, consistent, X, X_consistent = load()
    splits, assignments = split_manifest(y, groups, patients, rows)
    split_rows = [dict(**row, **{f"fold_{p}": assignments[(p, j)] for p in splits})
                  for j, row in enumerate(rows)]
    write_csv(OUT / "sample-identity-and-folds.csv", list(split_rows[0]), split_rows)

    raw_feature_rows = []
    for gene in genes:
        for cohort, data in frames.items():
            accs = data["by_gene"][gene]
            first = data["by_accession"][accs[0]][2]
            raw_feature_rows.append(dict(gene=gene, cohort=cohort, n_accessions=len(accs),
                                         first_accession=first, first_stable_accession=accs[0],
                                         all_stable_accessions="|".join(accs),
                                         old_first_accession_consistent=gene in consistent))
    write_csv(OUT / "gene-accession-audit.csv", list(raw_feature_rows[0]), raw_feature_rows)
    feature_rows = [dict(stable_accession=a, gene=frames[c]["by_accession"][a][1],
                         cohort=c, original_accession=frames[c]["by_accession"][a][2])
                    for a in accessions for c in COHORTS]
    write_csv(OUT / "primary-accession-map.csv", list(feature_rows[0]), feature_rows)

    records = []
    predict(X, y, rows, splits, "accession", records)
    predict(X_consistent, y, rows, splits, "consistent_gene_first_sensitivity", records)
    predict(X, y, rows, splits, "global_summaries", records)
    predict(X, y, rows, splits, "missingness", records)
    summary = summarize(records, rows)
    write_csv(OUT / "predictions.csv", list(records[0]), records)
    source_paths = [Path(cptac.__file__).resolve().parent / "data" /
                    f"umich-{c.lower()}" / "Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz"
                    for c in COHORTS]
    manifest = dict(status="candidate analysis, not a manuscript replacement",
                    feature_primary="shared stable ENSP accession, stripped of version; common gene annotation required",
                    sensitivity="old shared gene symbols with first stable accession consistent across cohorts",
                    tissue=".N adjacent normal and HNSCC .C cored normal are y=0; unsuffixed tumor y=1",
                    values="CPTAC UMICH reference-intensity-subtracted abundance; not raw total signal",
                    cohort_order=list(COHORTS), n_samples=len(rows), n_patients=len(set(patients)),
                    tumor=int(y.sum()), normal=int(len(y)-y.sum()), n_primary_features=len(accessions),
                    n_old_shared_gene_symbols=len(genes), n_consistent_first_gene_sensitivity=len(consistent),
                    fold_policy="StratifiedKFold random_state=28; true-patient GroupKFold 5; cohort LOGO",
                    model=XGB, software={k: version(k) for k in
                                        ("cptac", "xgboost", "numpy", "pandas", "scikit-learn")},
                    input_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
                    results=summary)
    (OUT / "candidate-results.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("saved", OUT, flush=True)


if __name__ == "__main__":
    main()
