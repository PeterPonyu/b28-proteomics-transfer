#!/usr/bin/env python3
"""Independent current-article checks, never a scientific submission clearance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import fitz
import numpy as np
from sklearn.metrics import roc_auc_score, brier_score_loss, log_loss

CAP = Path(__file__).resolve().parents[2]
MAN = CAP / "manuscript"
BASE = CAP / "analysis/revision_20260922/candidate-v1"
SIDE = BASE.parent / "candidate-v2-side-studies"
OUT = CAP / "output/verification"
RECORDED = CAP / "analysis/integration_20260923"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def current_build_inputs():
    """Every current TeX/table/bibliography/figure input, excluding historical packages."""
    paths = [MAN/"manuscript.tex", MAN/"SI.tex", MAN/"ref.bib"]
    paths += list((MAN/"sections").glob("*.tex")) + list((MAN/"tables").glob("*.tex"))
    paths += list((MAN/"figures").glob("b28_*.pdf"))
    return sorted(paths)


def current_generated_assets():
    return sorted([*MAN.glob("figures/b28_*.pdf"), *MAN.glob("figures/b28_*.png"),
                   *MAN.glob("figures/plotdata/*.csv"), *MAN.glob("tables/*.tex")])


def verify_build_identity(record):
    expected = {str(p.relative_to(CAP)) for p in current_build_inputs()}
    assert set(record["inputs_sha256"]) == expected, "Compiled input inventory changed; rebuild required"
    for name, digest in record["inputs_sha256"].items():
        assert sha(CAP/name) == digest, f"Compiled input changed; rebuild required: {name}"
    assert set(record["documents"]) == {"manuscript.pdf", "SI.pdf"}
    for name, doc in record["documents"].items():
        assert sha(MAN/name) == doc["sha256"], f"Compiled PDF changed; rebuild required: {name}"


def verify_layout_identity(record):
    expected = {str(p.relative_to(CAP)) for p in
                [MAN/"manuscript.pdf", MAN/"SI.pdf", *sorted(MAN.glob("figures/b28_*.pdf"))]}
    assert record["errors"] == [], "Automated layout defects are unresolved"
    assert record["rendered_in_this_run"] is True, "Layout record did not render current pages"
    rows = record["records"]
    assert len(rows) == len(expected) and {r["path"] for r in rows} == expected, \
        "Layout input inventory changed; render again"
    for item in rows:
        assert sha(CAP/item["path"]) == item["sha256"], f"Layout PDF changed: {item['path']}"
        with fitz.open(CAP/item["path"]) as pdf:
            assert len(item["pages"]) == len(item["rendered_files"]) == len(pdf), \
                f"Incomplete rendered-page inventory: {item['path']}"
        for raster in item["rendered_files"]:
            assert sha(CAP/raster["path"]) == raster["sha256"], f"Rendered page changed: {raster['path']}"


def metric(rows, field="p"):
    y = np.array([int(r["y"]) for r in rows])
    p = np.array([float(r[field]) for r in rows])
    bins = np.minimum(np.floor(p*10).astype(int), 9)
    ece = sum(np.mean(bins == b)*abs(np.mean(p[bins == b])-np.mean(y[bins == b]))
              for b in range(10) if np.any(bins == b))
    return {"auroc": float(roc_auc_score(y, p)), "ece_10": float(ece),
            "brier": float(brier_score_loss(y, p)), "log_loss": float(log_loss(y, p)),
            "accuracy": float(np.mean((p >= .5) == y))}


def verify_numerics():
    """Full independent numeric check. Ordinary PDF rebuilds reuse its pinned record."""
    identity = read(BASE/"sample-identity-and-folds.csv")
    assert len(identity) == len({(r["cohort"], r["sample_id"]) for r in identity}) == 1191
    assert len({r["patient_id"] for r in identity}) == 719
    assert Counter(r["tissue"] for r in identity) == {"tumor":689, "adjacent_normal":496, "cored_normal":6}
    for r in identity:
        sid = r["sample_id"]
        is_normal = sid.endswith(".N") or (r["cohort"] == "Hnscc" and sid.endswith(".C"))
        assert int(r["y"]) == int(not is_normal)
        assert r["patient_id"] == r["cohort"]+"|"+(sid[:-2] if is_normal else sid)
    crossing = {}
    for p in ("sample_random", "patient_grouped", "leave_one_cohort_out"):
        folds = defaultdict(set)
        for r in identity:
            folds[r["patient_id"]].add(r[f"fold_{p}"])
        crossing[p] = sum(len(s) > 1 for s in folds.values())
    assert crossing == {"sample_random":381, "patient_grouped":0, "leave_one_cohort_out":0}
    mapping = read(BASE/"primary-accession-map.csv")
    accs = defaultdict(list)
    for r in mapping:
        accs[r["stable_accession"]].append(r)
    assert len(accs) == 8906
    assert all(len(v)==6 and len({r["cohort"] for r in v})==6 and
               len({r["gene"] for r in v})==1 and len({r["original_accession"] for r in v})==1 for v in accs.values())
    groups = defaultdict(list)
    for r in read(BASE/"predictions.csv"):
        groups[r["feature_set"], r["protocol"]].append(r)
    results = json.loads((BASE/"candidate-results.json").read_text())["results"]
    cells = 0
    for (feature, protocol), rows in groups.items():
        assert len(rows) == 1191
        for k,v in metric(rows).items():
            assert abs(v-results[feature][protocol]["pooled"][k]) < 1e-12
            cells += 1
        for c in results[feature][protocol]["by_cohort"]:
            sub = [r for r in rows if r["cohort"] == c]
            for k,v in metric(sub).items():
                assert abs(v-results[feature][protocol]["by_cohort"][c][k]) < 1e-12
                cells += 1
    side_result = json.loads((SIDE/"results.json").read_text())
    side_pred = read(SIDE/"predictions.csv")
    for arm, data in side_result["results"].items():
        rows = [r for r in side_pred if r["arm"] == arm]
        assert len(rows) == 1191
        for c in (*data["by_cohort"], "pooled"):
            sub = rows if c == "pooled" else [r for r in rows if r["cohort"] == c]
            stored = data["pooled"] if c == "pooled" else data["by_cohort"][c]
            for k,v in metric(sub).items():
                assert abs(v-stored[k]) < 1e-12
                cells += 1
    assert side_result["train_median_l2"]["max_iter"] == 500
    assert all(r["converged"] and r["iterations"] < 500 for r in side_result["fit_details"] if r["arm"] == "train_median_l2")
    bound = json.loads((MAN/"integration-data.json").read_text())
    for field in ("source_sha256", "generated_sha256"):
        for name, expected in bound[field].items():
            assert sha(CAP/name) == expected, name
    transfer = read(MAN/"figures/plotdata/b28_transfer.csv")
    for row in transfer:
        c = row["cohort"]
        for column,key in (("auc_abundance","accession"), ("auc_missingness","missingness"), ("auc_trivial","global_summaries")):
            assert abs(float(row[column])-results[key]["leave_one_cohort_out"]["by_cohort"][c]["auroc"]) < 1e-12
    rec = read(BASE/"candidate-recalibration.csv")
    rec_summary = read(MAN/"figures/plotdata/b28_recal_summary.csv")
    for row in rec_summary:
        sub = rec if row["cohort"] == "pooled" else [r for r in rec if r["cohort"] == row["cohort"]]
        for col,field in (("ece_before","p_before"),("ece_after_temperature","p_temperature"),("ece_after_platt","p_platt")):
            assert abs(metric(sub,field)["ece_10"]-float(row[col])) < 1e-12
            cells += 1
    source_lookup = {(r["cohort"],r["sample_id"]):r for r in groups["accession","leave_one_cohort_out"]}
    ordered = [source_lookup[r["cohort"],r["sample_id"]] for r in identity]
    ranked = sorted(enumerate(ordered),key=lambda ir:(-abs(float(ir[1]["p"])-.5),ir[0]))
    for row in read(BASE/"candidate-riskcoverage.csv"):
        selected = [r for _,r in ranked[:int(row["n_retained"])]]
        assert abs(np.mean([(float(r["p"])>=.5)==int(r["y"]) for r in selected])-float(row["accuracy_retained"])) < 1e-12
        counts = Counter(r["cohort"] for r in selected)
        assert all(counts[c] == int(row[f"n_{c}"]) for c in results["accession"]["leave_one_cohort_out"]["by_cohort"])
    return dict(source_prediction_sets=len(groups), n_metric_cells=cells, true_patient_crossings=crossing)


def verify_numeric_identity(record):
    """Fail closed when a source/table used by the prior numeric check changes."""
    assert record["n_metric_cells"] == 301 and record["source_prediction_sets"] == 6
    if "prior_qa_path" in record:
        assert sha(CAP/record["prior_qa_path"]) == record["prior_qa_sha256"], "Prior numerical QA changed"
    for name, expected in record["verified_dependency_sha256"].items():
        assert sha(CAP/name) == expected, f"Numerical dependency changed; full verification required: {name}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse-numerics", action="store_true", help="reuse 301-cell verification only if every pinned dependency is unchanged")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    compilation_path = OUT/"compile-QA.json"
    if not compilation_path.exists():
        compilation_path = CAP/"provenance/public-compilation.json"
    if not compilation_path.exists():
        compilation_path = RECORDED/"compile-QA.json"
    compilation = json.loads(compilation_path.read_text())
    verify_build_identity(compilation)
    layout_path = OUT/"layout-QA.json"
    layout = json.loads(layout_path.read_text())
    verify_layout_identity(layout)
    recorded_numeric_path = RECORDED/"numerical-QA.json"
    # Archived dependency identity remains mandatory even on fresh recomputation.
    verify_numeric_identity(json.loads(recorded_numeric_path.read_text()))
    numeric_path = recorded_numeric_path
    if args.reuse_numerics:
        numeric = json.loads(numeric_path.read_text())
        verify_numeric_identity(numeric)
    else:
        numeric = verify_numerics()
        bound = json.loads((MAN/"integration-data.json").read_text())
        dependencies = dict(bound["source_sha256"], **bound["generated_sha256"])
        for p in (BASE/"primary-accession-map.csv", BASE/"gene-accession-audit.csv", SIDE/"predictions.csv"):
            dependencies[str(p.relative_to(CAP))] = sha(p)
        numeric.update(status="301 numerical metric cells independently recomputed",
                       verified_dependency_sha256=dependencies)
        numeric_path = OUT/"numerical-QA.json"
        numeric_path.write_text(json.dumps(numeric,indent=2)+"\n")
    bound = json.loads((MAN/"integration-data.json").read_text())
    for field in ("source_sha256", "generated_sha256"):
        for name, expected in bound[field].items():
            assert sha(CAP/name) == expected, name
    text = "\n".join(p.read_text() for p in [MAN/"manuscript.tex", *sorted((MAN/"sections").glob("*.tex"))])
    for term in ("0.946", "0.857--0.991", "73\\%", "27\\%", "116/66", "31.2\\%", "0.794", "0.831"):
        assert term not in text, f"Stale claim: {term}"
    for term in ("0.9343", "0.9555", "0.9942", "0.9944", "0.1100", "719", "8{,}906", "119/119"):
        assert term in text, term
    documents = {}
    for name in ("manuscript.pdf", "SI.pdf"):
        doc = fitz.open(MAN/name)
        words = "\n".join(p.get_text() for p in doc)
        assert "??" not in words
        assert len(doc) > (15 if name == "manuscript.pdf" else 4)
        if name == "manuscript.pdf":
            assert all(token in words for token in ("0.9343","0.9555","719","8,906"))
            assert "0.946" not in words and "73%" not in words
        assert not any(len(p.get_text().strip()) < 5 for p in doc)
        documents[name] = dict(pages=len(doc), sha256=sha(MAN/name))
    for p in (MAN/"figures").glob("b28_*.pdf"):
        doc = fitz.open(p)
        assert len(doc) == 1
        words = doc[0].get_text()
        assert not any(x in words for x in ("9675","9,675","0.946","one cancer type + one batch","no PDAC")), p
    # Original machine-specific logs are excluded. Check logs only after a local
    # explicit build; the preserved compile record binds the supplied PDFs.
    logs_checked = []
    for name in ("manuscript.log", "SI.log"):
        if (MAN/name).exists():
            log = (MAN/name).read_text(errors="replace")
            assert "undefined" not in log.lower() and "Overfull" not in log and "! " not in log, name
            logs_checked.append(name)
    ov_check = CAP/"provenance/historical/capsule-ov-source-check.json"
    ov_ledger = CAP/"provenance/historical/ov-summary-20260911.json"
    assert json.loads(ov_check.read_text())["ledger_sha256"] == sha(ov_ledger), \
        "Historical OV summary changed after source check"
    authored = [MAN/"manuscript.tex", MAN/"SI.tex", MAN/"ref.bib"]
    authored += sorted((MAN/"sections").glob("*.tex"))
    authored += sorted((MAN/"scripts").glob("*.py")) + sorted((MAN/"scripts").glob("*.sh"))
    authored += sorted((MAN/"scripts/R").glob("*.R"))
    authored += sorted((MAN/"scripts/tests").glob("*.py"))
    authored += [CAP/"reproducibility/verify_method_contract.py", CAP/"reproducibility/analysis/fold_identity.py",
                 CAP/"reproducibility/tests/test_b28_method_integrity.py"]
    report = dict(status="original article/SI integration checks passed; scientific HOLD", scientific_ready_attested=False,
                  source_prediction_sets=numeric["source_prediction_sets"], n_metric_cells=numeric["n_metric_cells"],
                  true_patient_crossings=numeric["true_patient_crossings"],
                  numerical_verification="reused with exact dependency hashes" if args.reuse_numerics else "recomputed",
                  numerical_record_sha256=sha(numeric_path),
                  compilation_record_sha256=sha(compilation_path),
                  layout_record_sha256=sha(layout_path),
                  verified_source_sha256={str(p.relative_to(CAP)):sha(p) for p in authored},
                  verified_asset_sha256={str(p.relative_to(CAP)):sha(p) for p in current_generated_assets()},
                  documents=documents, historical_submission_excluded=True,
                  upstream_verified_in_this_run=False, build_logs_checked=logs_checked,
                  ov_source_check=json.loads(ov_check.read_text()),
                  visual_review="Separate rendered-page inspection required; this verifier does not certify layout")
    (OUT/"integration-QA.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
