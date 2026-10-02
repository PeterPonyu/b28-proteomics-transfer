"""B28 candidate sensitivity protocols must preserve fold boundaries."""
import numpy as np
import pytest

from side_studies import fill_with_median, predict_arm


def test_only_deliberately_leaky_arm_sees_held_out_values():
    X = np.array([[np.nan, 1.], [2., 2.], [10., 100.]])
    tr, te = np.array([0, 1]), np.array([2])
    safe, _ = fill_with_median(X, tr, te, leak=False)
    leak, _ = fill_with_median(X, tr, te, leak=True)
    assert safe[0, 0] == 2.
    assert leak[0, 0] == 6.
    X[2, 0] = 1000.
    assert fill_with_median(X, tr, te, leak=False)[0][0, 0] == safe[0, 0]
    assert fill_with_median(X, tr, te, leak=True)[0][0, 0] != leak[0, 0]


def test_all_missing_train_column_gets_zero_not_held_out_value():
    X = np.array([[np.nan, 0.], [np.nan, 1.], [99., 3.]])
    tr, te = np.array([0, 1]), np.array([2])
    train, test = fill_with_median(X, tr, te, leak=False)
    assert list(train[:, 0]) == [0., 0.]
    assert test[0, 0] == 99.


def test_unknown_arm_fails_closed():
    with pytest.raises(ValueError, match="Unknown arm"):
        predict_arm("ambiguous", np.ones((4, 2)), np.array([0, 1, 0, 1]),
                    np.array([0, 1]), np.array([2, 3]))
