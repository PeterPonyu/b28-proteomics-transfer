"""Regressions for the corrected Article; never an endorsement of submission readiness."""
import csv
import importlib.util
import hashlib
import json
from pathlib import Path

import fitz
import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

CAP = Path(__file__).resolve().parents[3]
MAN = CAP / "manuscript"
BASE = CAP / "analysis/revision_20260922/candidate-v1"
SIDE = BASE.parent / "candidate-v2-side-studies"


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def source():
    return "\n".join(p.read_text() for p in [MAN/"manuscript.tex", *sorted((MAN/"sections").glob("*.tex"))])


def test_historical_submission_is_excluded_and_old_manifest_is_historical():
    assert not (MAN/'submission').exists()
    original = json.loads((CAP/'provenance/historical/capsule-release-current.json').read_text())
    assert original['historical_submission']['synchronized'] is False
    assert hashlib.sha256((MAN/'manuscript.pdf').read_bytes()).hexdigest() != original['historical_submission']['sha256']


def test_current_article_structure_and_population():
    text = source()
    assert text.count(r"\begin{figure}") == 8
    assert all(name in text for name in (r"\section{Introduction}", r"\section{Materials and Methods}",
                                          r"\section{Results}", r"\section{Discussion}"))
    assert "8{,}906" in text and "719" in text and "689" in text and "502" in text
    assert "\\texttt{.C}" in text and "470" in text and "381" in text
    for forbidden in ("73\\%", "27\\%", "0.946", "0.857--0.991", "116/66", "31.2\\%"):
        assert forbidden not in text


def test_macro_is_cohort_mean_not_pooled_auc():
    rows = [r for r in read(BASE/"predictions.csv") if r["feature_set"]=="accession" and r["protocol"]=="leave_one_cohort_out"]
    cohort_auc=[]
    for c in sorted({r["cohort"] for r in rows}):
        sub=[r for r in rows if r["cohort"]==c]
        cohort_auc.append(roc_auc_score([int(r["y"]) for r in sub],[float(r["p"]) for r in sub]))
    macro=float(np.mean(cohort_auc))
    pooled=roc_auc_score([int(r["y"]) for r in rows],[float(r["p"]) for r in rows])
    assert macro == pytest.approx(.9343426660543135,abs=1e-12)
    assert pooled == pytest.approx(.9555305627995999,abs=1e-12)
    assert abs(macro-pooled) > .02
    assert "0.9343" in source() and "0.9555" in source()


def test_l2_executed_iteration_limit_and_convergence():
    side=json.loads((SIDE/"results.json").read_text())
    assert side["train_median_l2"]["max_iter"] == 500
    fits=[r for r in side["fit_details"] if r["arm"]=="train_median_l2"]
    assert len(fits)==6 and all(r["converged"] and r["iterations"]<500 for r in fits)
    assert r"\texttt{max\_iter=500}" in source()


def test_corrected_hnscc_counts_and_coverage_sources():
    h=next(r for r in read(MAN/"figures/plotdata/b28_pair_counts.csv") if r["cohort"]=="Hnscc")
    assert (int(h["n_tumor"]),int(h["n_normal"]),int(h["n_patients"]),int(h["n_paired_patients"])) == (110,72,112,68)
    low=read(BASE/"candidate-riskcoverage.csv")[0]
    assert [int(low[f"n_{c}"]) for c in ("Ccrcc","Luad","Ucec","Hnscc","Lscc","Pdac")] == [6,52,0,4,55,2]
    assert int(low["n_retained"])==119 and float(low["accuracy_retained"])==1


def test_pdf_text_and_references_match_new_source():
    for name in ("manuscript", "SI"):
        doc=fitz.open(MAN/f"{name}.pdf")
        text="\n".join(p.get_text() for p in doc)
        assert "??" not in text
        assert "719" in text and "8,906" in text
        assert "0.946" not in text and "73%" not in text
        if (MAN/f"{name}.log").exists():
            log=(MAN/f"{name}.log").read_text()
            assert "Overfull" not in log and "undefined" not in log.lower()


def test_current_release_never_claims_scientific_or_submission_sync():
    release=json.loads((CAP/"PUBLIC-MANIFEST.json").read_text())
    assert release["scientific_ready_attested"] is False
    assert release["historical_submission_included"] is False
    assert "HOLD" in release["status"]
    assert {r["path"] for r in release["documents"]} == {"manuscript/manuscript.pdf","manuscript/SI.pdf"}
    assert release["schema_version"] == "b28-public-current-1"


def test_offline_build_never_invokes_scientific_fit_or_historical_exporter():
    script=(MAN/"scripts/rebuild_current.sh").read_text()
    assert "export_revision.py" in script and "verify_integration.py" in script
    assert "python analysis/revision_20260922/reanalyse.py" not in script
    assert "python analysis/revision_20260922/side_studies.py" not in script
    assert "compile_current.py" in script and "check_layout.py --render" in script
    assert '"SI"' in (MAN/"scripts/compile_current.py").read_text()
    assert "--reuse-numerics" in script and "verify_method_contract.py" in script
    assert '--build' in script and 'release_integration.py' not in script
    assert "results-snapshot/B28_real" not in (MAN/"scripts/make_study_design.py").read_text()


def test_every_side_study_has_explicit_current_disposition():
    text=(MAN/"SI.tex").read_text()
    for term in ("OV plex-LOGO", "label-permutation", "ComBat", "PC-removal", "Residual", "injected-shift",
                 "Two-way variance", "BCM/UMICH", "BRCA", "COAD", "GBM"):
        assert term in text
    assert "not promoted" in text
    assert "summary-level historical results" in text
    assert "not a new prediction-level replay" in " ".join(text.split())


def test_reused_metric_verification_rejects_modified_dependency(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("b28_integration_check", MAN/"scripts/verify_integration.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "CAP", tmp_path)
    data = tmp_path/"sample.csv"
    data.write_text("sample,value\ns1,1\n")
    record = dict(n_metric_cells=301, source_prediction_sets=6,
                  verified_dependency_sha256={"sample.csv":hashlib.sha256(data.read_bytes()).hexdigest()})
    module.verify_numeric_identity(record)
    data.write_text("sample,value\ns1,2\n")
    with pytest.raises(AssertionError, match="Numerical dependency changed"):
        module.verify_numeric_identity(record)


def test_current_contract_rejects_historical_schema():
    contract=json.loads((CAP/"reproducibility/method-contract.json").read_text())
    assert contract["schema_version"] == "b28-public-current-1"
    paths={r["path"] for r in contract["locked_files"]}
    assert "analysis/revision_20260922/candidate-v1/predictions.csv" in paths
    assert "reproducibility/splits/patient-grouped-corrected-20260922.csv" in paths
    assert not any(p.startswith("results-snapshot/B28_real/") for p in paths)


def integration_module():
    spec = importlib.util.spec_from_file_location("b28_artifact_check", MAN/"scripts/verify_integration.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_compilation_guard_rejects_source_changed_without_rebuild(tmp_path, monkeypatch):
    module = integration_module()
    monkeypatch.setattr(module, "CAP", tmp_path)
    monkeypatch.setattr(module, "MAN", tmp_path)
    source = tmp_path/"manuscript.tex"
    source.write_text("Current corrected manuscript")
    for name in ("manuscript.pdf", "SI.pdf"):
        (tmp_path/name).write_bytes(b"compiled PDF bytes for hash test")
    monkeypatch.setattr(module, "current_build_inputs", lambda: [source])
    record = dict(inputs_sha256={"manuscript.tex":module.sha(source)},
                  documents={name:dict(sha256=module.sha(tmp_path/name)) for name in ("manuscript.pdf", "SI.pdf")})
    module.verify_build_identity(record)
    source.write_text("Changed text not yet compiled")
    with pytest.raises(AssertionError, match="Compiled input changed"):
        module.verify_build_identity(record)


def test_compilation_guard_rejects_replaced_pdf(tmp_path, monkeypatch):
    module = integration_module()
    monkeypatch.setattr(module, "CAP", tmp_path)
    monkeypatch.setattr(module, "MAN", tmp_path)
    monkeypatch.setattr(module, "current_build_inputs", lambda: [])
    for name in ("manuscript.pdf", "SI.pdf"):
        (tmp_path/name).write_bytes(b"current compiled bytes")
    record = dict(inputs_sha256={}, documents={name:dict(sha256=module.sha(tmp_path/name))
                                              for name in ("manuscript.pdf", "SI.pdf")})
    module.verify_build_identity(record)
    (tmp_path/"SI.pdf").write_bytes(b"superseded SI")
    with pytest.raises(AssertionError, match="Compiled PDF changed"):
        module.verify_build_identity(record)


def test_layout_guard_rejects_modified_render(tmp_path, monkeypatch):
    module = integration_module()
    monkeypatch.setattr(module, "CAP", tmp_path)
    monkeypatch.setattr(module, "MAN", tmp_path)
    records = []
    for name in ("manuscript", "SI"):
        pdf = fitz.open()
        pdf.new_page().insert_text((30,30), "Render identity test")
        pdf.save(tmp_path/f"{name}.pdf")
        pdf.close()
        raster = tmp_path/f"{name}.png"
        raster.write_bytes(b"raster bytes for hash test")
        records.append(dict(path=f"{name}.pdf", sha256=module.sha(tmp_path/f"{name}.pdf"),
                            pages=[dict(page=1)], rendered_files=[dict(path=raster.name,sha256=module.sha(raster))]))
    record = dict(errors=[], records=records, rendered_in_this_run=True)
    module.verify_layout_identity(record)
    (tmp_path/"SI.png").write_bytes(b"old raster substituted")
    with pytest.raises(AssertionError, match="Rendered page changed"):
        module.verify_layout_identity(record)


def test_release_hashes_cover_current_files_and_do_not_inherit_visual_pass():
    release = json.loads((CAP/"PUBLIC-MANIFEST.json").read_text())
    for item in release['locked_files'] + release['control_files']:
        assert hashlib.sha256((CAP/item['path']).read_bytes()).hexdigest() == item['sha256'], item['path']
    assert release['visual_review']['scientific_review_completed'] is False
    assert release['scientific_ready_attested'] is False
    for item in release['documents']:
        with fitz.open(CAP/item['path']) as pdf:
            assert len(pdf) == item['pages']


def test_method_contract_actually_rejects_old_schema_and_changed_evidence(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("b28_method_check", CAP/"reproducibility/verify_method_contract.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with pytest.raises(SystemExit, match="Historical method contract"):
        module.validate(dict(schema_version="historical"))
    monkeypatch.setattr(module, "CAP", tmp_path)
    paths = ("analysis/revision_20260922/candidate-v1/predictions.csv",
             "analysis/revision_20260922/candidate-v2-side-studies/predictions.csv",
             "analysis/integration_20260923/numerical-QA.json")
    items = []
    for name in paths:
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"unchanged corrected evidence")
        items.append(dict(path=name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    record = dict(paper="B28",schema_version="b28-public-current-1", scientific_ready_attested=False,
                  status="SCIENTIFIC_HOLD",locked_files=items)
    module.validate(record, evidence_only=True)
    (tmp_path/paths[0]).write_bytes(b"old or modified predictions")
    with pytest.raises(SystemExit, match="SOURCE IDENTITY CHANGED"):
        module.validate(record, evidence_only=True)
