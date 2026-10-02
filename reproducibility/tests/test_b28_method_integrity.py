"""The corrected patient partition follows true identity, never row position."""
import ast
import csv
from pathlib import Path
import sys
import warnings
import numpy as np
import pytest

REPRO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO/'analysis'))
from fold_identity import replay_patient_splits


@pytest.fixture
def rows():
    with (REPRO/'splits/patient-grouped-corrected-20260922.csv').open() as f:
        return list(csv.DictReader(f))


def replay(locked, presented=None):
    r = locked if presented is None else presented
    return replay_patient_splits(locked, [v['cohort'] for v in r], [v['sample_id'] for v in r],
                                 [v['patient_key'] for v in r], [int(v['label']) for v in r])


def test_real_partition_is_exhaustive_patient_disjoint_and_order_invariant(rows):
    assert len(rows) == 1191 and len({r['patient_key'] for r in rows}) == 719
    a = replay(rows)
    shuffled = [rows[i] for i in np.random.RandomState(7).permutation(len(rows))]
    b = replay(rows, shuffled)
    assert sorted(np.concatenate([test for _, test in a])) == list(range(len(rows)))
    for (tr, te), (_, shuffled_te) in zip(a, b):
        assert {rows[i]['patient_key'] for i in tr}.isdisjoint({rows[i]['patient_key'] for i in te})
        identity = lambda r: (r['cohort'], r['sample_id'])
        assert {identity(rows[i]) for i in te} == {identity(shuffled[i]) for i in shuffled_te}


def test_duplicate_or_missing_sample_is_rejected(rows):
    with pytest.raises(ValueError, match='one-to-one'):
        replay(rows, rows[:-1]+[rows[0]])
    with pytest.raises(ValueError, match='inventory'):
        replay(rows, rows[:-1])


@pytest.mark.parametrize('change', [{'label': '0'}, {'patient_key': 'new-patient'}])
def test_changed_label_or_patient_is_rejected(rows, change):
    with pytest.raises(ValueError, match='identity changed'):
        replay(rows, [dict(rows[0], **change)]+rows[1:])


def test_locked_table_cannot_split_one_patient(rows):
    paired = next(r['patient_key'] for r in rows if sum(x['patient_key']==r['patient_key'] for x in rows)>1)
    idx = next(i for i, r in enumerate(rows) if r['patient_key']==paired)
    rows[idx] = dict(rows[idx], fold=str(int(rows[idx]['fold']) % 5 + 1))
    with pytest.raises(ValueError, match='patient crosses'):
        replay(rows)


def test_hnscc_cored_normals_share_true_patient_and_fold(rows):
    cored = [r for r in rows if r['cohort']=='Hnscc' and r['sample_id'].endswith('.C')]
    assert len(cored) == 6
    for row in cored:
        assert row['label'] == '0'
        assert row['patient_key'] == 'Hnscc|'+row['sample_id'][:-2]
        assert len({r['fold'] for r in rows if r['patient_key']==row['patient_key']}) == 1


def test_held_out_values_cannot_change_training_imputation():
    # Execute the corrected side-study preprocessor, not the archived matrix pipeline.
    tree = ast.parse((REPRO.parent/'analysis/revision_20260922/side_studies.py').read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name=='fill_with_median')
    namespace = {'np': np, 'warnings': warnings}
    exec(compile(ast.Module(body=[function], type_ignores=[]), '<impute_train>', 'exec'), namespace)
    train = np.array([[1., 2.], [3., np.nan]])
    held = np.array([[99999., np.nan], [np.nan, -99999.]])
    X = np.vstack((train, held))
    out_train, out_test = namespace['fill_with_median'](X, np.array([0,1]), np.array([2,3]), leak=False)
    np.testing.assert_array_equal(out_train, [[1, 2], [3, 2]])
    np.testing.assert_array_equal(out_test, [[99999, 2], [2, -99999]])
    changed = X.copy()
    changed[2:] = [[-1e8, -1e8], [1e8, 1e8]]
    changed_train, _ = namespace['fill_with_median'](changed, np.array([0,1]), np.array([2,3]), leak=False)
    np.testing.assert_array_equal(out_train, changed_train)
