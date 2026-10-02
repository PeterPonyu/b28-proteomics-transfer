"""Fail-closed public included-file closure, never scientific clearance."""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json

CAP = Path(__file__).resolve().parents[1]
SCHEMA = 'b28-public-current-1'


def validate(record, evidence_only=False):
    if record.get('schema_version') != SCHEMA:
        raise SystemExit('Historical method contract cannot validate the public export')
    if record.get('scientific_ready_attested') is not False or 'HOLD' not in record.get('status', ''):
        raise SystemExit('This build cannot certify scientific readiness')
    items = record['locked_files']
    names = [item['path'] for item in items]
    if not items or len(set(names)) != len(names):
        raise SystemExit('Empty or duplicate public dependency inventory')
    required = {'analysis/revision_20260922/candidate-v1/predictions.csv',
                'analysis/revision_20260922/candidate-v2-side-studies/predictions.csv',
                'analysis/integration_20260923/numerical-QA.json'}
    if not required.issubset(names):
        raise SystemExit('Corrected evidence is absent from the public contract')
    if evidence_only:
        # Source edits may be rebuilt; they are bound again only after actual compilation
        # and integration QA. Scientific evidence and its prior QA cannot be refreshed here.
        items = [item for item in items if item['path'].startswith('analysis/')]
    fail = []
    for item in items:
        name = PurePosixPath(item['path'])
        path = CAP/name
        if (name.is_absolute() or '..' in name.parts or path.is_symlink()
                or not path.resolve().is_relative_to(CAP.resolve()) or not path.is_file()
                or hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']):
            fail.append(item['path'])
    if fail:
        raise SystemExit('SOURCE IDENTITY CHANGED; review protocol before rebuilding: '+', '.join(fail))
    if not evidence_only:
        current = public_files()
        if current != set(names):
            raise SystemExit('PUBLIC INVENTORY CHANGED: '+', '.join(sorted(current ^ set(names))))
    print(f"B28: {len(items)} included files match; SCIENTIFIC_HOLD; upstream matrix not checked")


def public_files():
    controls = {'PUBLIC-MANIFEST.json', 'SHA256SUMS.txt', 'reproducibility/method-contract.json'}
    transient_dirs = {'output', '.git', '__pycache__', '.pytest_cache'}
    transient_suffixes = {'.aux', '.bbl', '.blg', '.fdb_latexmk', '.fls', '.log', '.out', '.pyc'}
    return {p.relative_to(CAP).as_posix() for p in CAP.rglob('*') if p.is_file()
            and p.relative_to(CAP).as_posix() not in controls
            and not transient_dirs.intersection(p.relative_to(CAP).parts)
            and p.suffix not in transient_suffixes}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence-only', action='store_true',
                        help='preflight unchanged scientific evidence before rebuilding source edits')
    validate(json.loads((CAP/'reproducibility/method-contract.json').read_text()),
             parser.parse_args().evidence_only)
