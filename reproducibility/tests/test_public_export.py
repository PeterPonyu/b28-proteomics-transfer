"""Public path/closure and offline-boundary regressions; no fitted models."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
REQUIRED = ('analysis/revision_20260922/candidate-v1/predictions.csv',
            'analysis/revision_20260922/candidate-v2-side-studies/predictions.csv',
            'analysis/integration_20260923/numerical-QA.json')


def module():
    spec = importlib.util.spec_from_file_location('public_guard', ROOT/'reproducibility/verify_method_contract.py')
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def minimal(tmp_path):
    entries = []
    for name in REQUIRED:
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'unchanged scientific evidence')
        entries.append(dict(path=name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return dict(paper='B28', schema_version='b28-public-current-1', status='SCIENTIFIC_HOLD',
                scientific_ready_attested=False, locked_files=entries)


def test_default_guard_rejects_changed_and_missing_evidence(tmp_path, monkeypatch):
    guard = module()
    monkeypatch.setattr(guard, 'CAP', tmp_path)
    record = minimal(tmp_path)
    guard.validate(record)
    path = tmp_path/REQUIRED[0]
    path.write_bytes(b'changed predictions')
    with pytest.raises(SystemExit, match='SOURCE IDENTITY CHANGED'):
        guard.validate(record)
    path.unlink()
    with pytest.raises(SystemExit, match='SOURCE IDENTITY CHANGED'):
        guard.validate(record)


def test_guard_rejects_traversal_symlink_and_unexpected_file(tmp_path, monkeypatch):
    guard = module()
    monkeypatch.setattr(guard, 'CAP', tmp_path)
    record = minimal(tmp_path)
    record['locked_files'].append(dict(path='../escape.csv', sha256='0'*64))
    with pytest.raises(SystemExit, match='SOURCE IDENTITY CHANGED'):
        guard.validate(record)
    record['locked_files'].pop()
    path = tmp_path/REQUIRED[0]
    actual = tmp_path/'link-target'
    actual.write_bytes(path.read_bytes())
    path.unlink()
    path.symlink_to(actual)
    with pytest.raises(SystemExit, match='SOURCE IDENTITY CHANGED'):
        guard.validate(record)
    path.unlink()
    path.write_bytes(actual.read_bytes())
    with pytest.raises(SystemExit, match='PUBLIC INVENTORY CHANGED'):
        guard.validate(record)


def test_default_actual_tree_passes_and_output_is_not_payload():
    guard = module()
    record = json.loads((ROOT/'reproducibility/method-contract.json').read_text())
    guard.validate(record)
    assert not any(p.startswith('output/') for p in guard.public_files())


def test_current_code_function_import_does_not_import_cptac():
    code = "import sys; sys.path.insert(0, sys.argv[1]); import reanalyse; assert 'cptac' not in sys.modules; assert reanalyse.stable_accession('ENSP00000432830.1') == 'ENSP00000432830'"
    subprocess.run([sys.executable, '-B', '-c', code, str(ROOT/'analysis/revision_20260922')], check=True)


def test_missing_optional_upstream_fails_instead_of_downloading(tmp_path):
    for name in ('verify_candidate.py', 'verify_side_studies.py'):
        result = subprocess.run([sys.executable, '-B', str(ROOT/'analysis/revision_20260922'/name),
                                 '--upstream-root', str(tmp_path)], capture_output=True, text=True)
        assert result.returncode != 0
        assert 'Upstream source absent or changed' in result.stderr


def test_payload_has_no_private_submission_backups_or_raw_proteomics():
    names = module().public_files()
    forbidden = {'submission', 'backups', '.omc', 'original', 'results-snapshot'}
    assert all(not forbidden.intersection(Path(n).parts) for n in names)
    assert not any('cover_letter' in n.lower() or n.endswith(('.bak','.tar.gz','.zip','.tsv.gz')) for n in names)
    assert not (ROOT/'SOURCE_PIN.md').exists()
    assert not (ROOT/'docs/index.html').exists()


def test_upstream_provenance_does_not_invent_pdc_accessions():
    data = json.loads((ROOT/'provenance/data-sources.json').read_text())
    assert len(data['cohorts']) == 6
    assert all(c['pdc_study_accession'] is None and c['source_file_in_public_release'] is False
               for c in data['cohorts'])


def test_version_identifier_does_not_self_attest_remote_publication():
    text = (ROOT/'README.md').read_text()
    assert '10.5281/zenodo.23101188' in text
    assert 'identifier alone is not proof of publication' in text
    assert 'https://github.com/PeterPonyu/b28-proteomics-transfer' in text
    public = json.loads((ROOT/'PUBLIC-MANIFEST.json').read_text())
    assert public['zenodo']['version_doi'] == '10.5281/zenodo.23101188'
    assert public['zenodo']['public_publication_verified'] is None
    assert public['scientific_ready_attested'] is False
