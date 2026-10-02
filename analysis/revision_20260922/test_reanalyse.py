"""Identity and estimand regressions for the candidate B28 rerun."""
import numpy as np
import pytest

from reanalyse import ece, split_manifest, stable_accession, tissue_patient


def test_known_normal_and_true_patient_identity():
    assert tissue_patient("Hnscc", "C3L-00994.C") == (
        "cored_normal", 0, "Hnscc|C3L-00994")
    assert tissue_patient("Hnscc", "C3L-00994.N") == (
        "adjacent_normal", 0, "Hnscc|C3L-00994")
    assert tissue_patient("Hnscc", "C3L-00994") == (
        "tumor", 1, "Hnscc|C3L-00994")
    with pytest.raises(ValueError, match="Unknown sample status"):
        tissue_patient("Hnscc", "C3L-00994.X")
    with pytest.raises(ValueError, match="Unknown sample status"):
        tissue_patient("Lscc", "C3L-00994.C")


def test_stable_accession_fail_closed():
    assert stable_accession("ENSP00000432830.1") == "ENSP00000432830"
    with pytest.raises(ValueError, match="Unsupported protein accession"):
        stable_accession("ENSP00000432830|ENSP00000377078")


def test_split_keeps_triplet_together_and_exhausts_rows():
    rows = [dict(cohort="A", sample_id=f"A{i}") for i in range(12)] + [
        dict(cohort="B", sample_id=f"B{i}") for i in range(12)]
    patients = np.array(["A|p0"]*3 + [f"A|p{i}" for i in range(1, 10)] +
                        [f"B|p{i}" for i in range(12)])
    y = np.array([0, 0, 1] + [0, 1]*4 + [0] + [0, 1]*6)
    groups = np.array([r["cohort"] for r in rows])
    splits, assign = split_manifest(y, groups, patients, rows)
    assert len(assign) == 3*len(rows)
    assert len({assign[("patient_grouped", i)] for i in range(3)}) == 1
    for protocol in ("patient_grouped", "leave_one_cohort_out"):
        for tr, te in splits[protocol]:
            assert set(patients[tr]).isdisjoint(set(patients[te]))


def test_ece_constant_prevalence_is_zero():
    y = np.array([0, 0, 1, 1])
    p = np.full(4, .5)
    assert ece(y, p) == 0.
