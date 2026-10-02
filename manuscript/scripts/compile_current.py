#!/usr/bin/env python3
"""Compile the original Article and SI and bind the successful build to its inputs."""
from __future__ import annotations

import json
import os
import subprocess

import fitz

from verify_integration import CAP, MAN, OUT, current_build_inputs, sha


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = {str(p.relative_to(CAP)): sha(p) for p in current_build_inputs()}
    environment = dict(os.environ, SOURCE_DATE_EPOCH="0", FORCE_SOURCE_DATE="1", TZ="UTC")
    for name in ("manuscript", "SI"):
        subprocess.run(["latexmk", "-g", "-pdf", "-interaction=nonstopmode",
                        "-halt-on-error", f"{name}.tex"], cwd=MAN, env=environment, check=True)
    assert inputs == {str(p.relative_to(CAP)): sha(p) for p in current_build_inputs()}, \
        "Build inputs changed during compilation; no build attestation written"
    documents = {}
    for name in ("manuscript.pdf", "SI.pdf"):
        with fitz.open(MAN/name) as pdf:
            documents[name] = dict(pages=len(pdf), sha256=sha(MAN/name))
    report = dict(status="Article and SI compiled successfully by latexmk",
                  scientific_ready_attested=False, inputs_sha256=inputs, documents=documents)
    (OUT/"compile-QA.json").write_text(json.dumps(report, indent=2)+"\n")
    print(f"Compilation bound {len(inputs)} source/figure inputs to Article and SI")


if __name__ == "__main__":
    main()
