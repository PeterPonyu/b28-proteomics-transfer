"""Explicit reviewed public inventory writer; not part of default reproduction.

Never rebind scientific outputs. A source/PDF change must have a matching actual
successful compilation record before the writer can record its new file hashes.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

from verify_method_contract import CAP, SCHEMA, public_files


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entry(name):
    path = CAP/name
    return dict(path=name, sha256=sha(path), bytes=path.stat().st_size)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-reviewed-inventory', action='store_true')
    args = parser.parse_args()
    if not args.write_reviewed_inventory:
        raise SystemExit('Explicit --write-reviewed-inventory required; never called by reproduce scripts')
    ledger = json.loads((CAP/'provenance/export-source-ledger.json').read_text())
    scientific = [r for r in ledger['files'] if r['path'].startswith(('analysis/','reproducibility/splits/'))
                  and r['path'].endswith(('.csv','.json'))]
    for row in scientific:
        if sha(CAP/row['path']) != row['source_sha256']:
            raise SystemExit('Frozen scientific output changed; cannot rebind: '+row['path'])
    spec = importlib.util.spec_from_file_location('b28_inventory_integration', CAP/'manuscript/scripts/verify_integration.py')
    check = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(check)
    check.verify_numeric_identity(json.loads((check.RECORDED/'numerical-QA.json').read_text()))
    compiled = check.OUT/'compile-QA.json'
    if not compiled.exists():
        compiled = CAP/'provenance/public-compilation.json'
    if not compiled.exists():
        compiled = check.RECORDED/'compile-QA.json'
    compile_record = json.loads(compiled.read_text())
    check.verify_build_identity(compile_record)
    if compiled.is_relative_to(check.OUT):
        # A new legitimate build record must join the shipped file closure,
        # not remain an unbundled dependency under disposable output/.
        public_compiled = CAP/'provenance/public-compilation.json'
        public_compiled.write_text(compiled.read_text())
        compiled = public_compiled
    names = sorted(public_files())
    if any((CAP/n).is_symlink() for n in names):
        raise SystemExit('Public payload may not contain symlinks')
    modified = [dict(path=r['path'], source_sha256=r['source_sha256'], public_sha256=sha(CAP/r['path']))
                for r in ledger['files'] if sha(CAP/r['path']) != r['source_sha256']]
    files = [entry(n) for n in names]
    contract = dict(paper='B28', schema_version=SCHEMA,
        status='PUBLIC_RESEARCH_EXPORT; SCIENTIFIC_HOLD', scientific_ready_attested=False,
        no_new_fit_scores_or_inference=True, locked_files=files,
        closure='Actual included files; original capsule manifests remain historical',
        upstream_raw_data_included=False)
    contract_path = CAP/'reproducibility/method-contract.json'
    contract_path.write_text(json.dumps(contract, indent=2)+'\n')
    docs = [dict(**entry('manuscript/'+name), pages=value['pages'])
            for name,value in compile_record['documents'].items()]
    manifest = dict(**contract, date='2026-10-02',
        repository='https://github.com/PeterPonyu/b28-proteomics-transfer',
        zenodo=dict(version_doi='10.5281/zenodo.23101188',
                    concept_doi='10.5281/zenodo.21870023',
                    tag='author-review-2026-10-02',
                    status='VERSION_IDENTIFIER_ASSIGNED; PUBLICATION_VERIFICATION_SEPARATE',
                    public_publication_verified=None,
                    original_concept_preserved=True),
        payload_files=len(files), payload_bytes=sum(r['bytes'] for r in files),
        control_files=[entry('reproducibility/method-contract.json')],
        documents=docs, historical_submission_included=False,
        numerical_dependency_count=39, numeric_metric_cells=301,
        source_ledger='provenance/export-source-ledger.json', modified_source_copies=modified,
        compilation_record=dict(path=compiled.relative_to(CAP).as_posix(), sha256=sha(compiled),
            boundary='Actual preserved or explicitly rebuilt compilation input identity; no original private dependency closure claim'),
        visual_review=dict(scientific_review_completed=False,
            status='Fresh export visual review recorded separately; historical visual PASS not inherited'),
        historical_records='provenance/historical/README.md')
    manifest_path = CAP/'PUBLIC-MANIFEST.json'
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    checksum_names = names+['reproducibility/method-contract.json','PUBLIC-MANIFEST.json']
    (CAP/'SHA256SUMS.txt').write_text(''.join(f'{sha(CAP/n)}  {n}\n' for n in sorted(checksum_names)))
    print(f'Public inventory: {len(files)} payload files; {len(scientific)} immutable scientific inputs; '
          f'{len(modified)} bounded source adaptations; HOLD')


if __name__=='__main__':
    main()
