"""Replay the audited cohort/sample split without silently changing patients."""
import numpy as np


def replay_patient_splits(rows, groups, sample_ids, patients, labels):
    n = len(labels)
    if not all(len(v) == n for v in (rows, groups, sample_ids, patients)):
        raise ValueError("locked split sample inventory changed")
    index = {(r['cohort'], r['sample_id']): r for r in rows}
    identities = [(str(c), str(s)) for c, s in zip(groups, sample_ids)]
    if len(index) != n or len(set(identities)) != n or set(index) != set(identities):
        raise ValueError("sample identities must match the locked inventory one-to-one")
    folds, patient_folds = [], {}
    for key, patient, label in zip(identities, patients, labels):
        row = index[key]
        if row['patient_key'] != patient or int(row['label']) != label:
            raise ValueError("patient/label identity changed")
        fold = int(row['fold'])
        if patient in patient_folds and patient_folds[patient] != fold:
            raise ValueError("one patient crosses fit/evaluation folds")
        patient_folds[patient] = fold
        folds.append(fold)
    folds = np.asarray(folds)
    if set(folds) != set(range(1, 6)):
        raise ValueError("expected all five locked folds")
    return [(np.flatnonzero(folds != i), np.flatnonzero(folds == i)) for i in range(1, 6)]
