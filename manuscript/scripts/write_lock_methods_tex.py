#!/usr/bin/env python3
"""Validate current executed method settings; never rewrite historical submission packs."""
import json
from pathlib import Path

CAP = Path(__file__).resolve().parents[2]
side = json.loads((CAP/"analysis/revision_20260922/candidate-v2-side-studies/results.json").read_text())
method = side["train_median_l2"]
assert method["max_iter"] == 500 and method["C"] == 1.0 and method["solver"] == "lbfgs"
assert all(r["converged"] for r in side["fit_details"] if r["arm"] == "train_median_l2")
print("Current methods: manuscript/sections/methods-current.tex and manuscript/SI.tex")
print("L2: train-fold median/scaling, C=1, lbfgs, max_iter=500; six converged fits")
print("Historical default-article methods-lock.tex is archived and not regenerated.")
