"""Shared atomic JSON, content hashes, and covariance provenance for swordfish runs."""
import hashlib
import json
import os
from pathlib import Path
import tempfile


def atomic_json(path, value):
    """Flush JSON to a sibling temporary file, then atomically replace one JSON artifact.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.
    value
        Value to serialize or include in the generated fixture record.

    Notes
    -----
    All hashes are content based. Covariance identity covers the complete native fit and relevant source/configuration, so touching timestamps cannot make a stale covariance valid. JSON publication uses a sibling temporary file and replacement; this is not a distributed-filesystem transaction.
    """
    # Keep the temporary sibling on the destination filesystem. Replace only
    # after serialization, flush and fsync succeed; failed JSON leaves the old
    # artifact intact. The caller coordinates larger multi-file transactions.
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix='.'+path.name, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def file_hash(path):
    """Return streaming SHA-256 for a file without loading large reports into memory.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.

    Returns
    -------
    Lowercase SHA-256 hex digest; file bytes are streamed in bounded chunks.

    Notes
    -----
    All hashes are content based. Covariance identity covers the complete native fit and relevant source/configuration, so touching timestamps cannot make a stale covariance valid. JSON publication uses a sibling temporary file and replacement; this is not a distributed-filesystem transaction.
    """
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def source_hashes(root, repo=None):
    """Fingerprint raw inputs, configuration, and case/framework source affecting a run.

    Parameters
    ----------
    root
        Example directory used as the reference point for input/output paths.
    repo
        Optional repository root; inferred from the supported example layout when omitted.

    Returns
    -------
    Mapping from repository-relative paths to content digests.

    Notes
    -----
    All hashes are content based. Covariance identity covers the complete native fit and relevant source/configuration, so touching timestamps cannot make a stale covariance valid. JSON publication uses a sibling temporary file and replacement; this is not a distributed-filesystem transaction.
    """
    root = Path(root).resolve()
    repo = Path(repo).resolve() if repo is not None else root.parents[2]
    paths = [root/'Makefile', root/'reference/source_manifest.json', repo/'Makefile']
    paths.extend((repo/'tests').glob('*.cpp'))
    paths.extend((root/'reference').glob('published*'))
    for directory in ['include', 'src', 'scripts', 'tests']:
        paths.extend(p for p in (root/directory).rglob('*') if p.suffix in {'.hpp', '.cpp', '.py', '.R'})
    paths.extend(p for p in (root/'data').rglob('*') if p.is_file())
    for directory in ['core', 'model', 'math', 'include', 'external/eigen', 'external/LBFGSpp']:
        paths.extend(p for p in (repo/directory).rglob('*') if p.is_file())
    result = {}
    for p in sorted(set(paths)):
        if '__pycache__' not in p.parts and p.name != '.DS_Store':
            result[str(p.relative_to(repo))] = file_hash(p)
    return result


def covariance_inputs(root, repo=None):
    """Bind curvature to all source/input hashes and the complete fitted native vector.

    Parameters
    ----------
    root
        Example directory used as the reference point for input/output paths.
    repo
        Optional repository root; inferred from the supported example layout when omitted.

    Returns
    -------
    Identity dictionary binding relevant source and the complete fitted native vector.

    Notes
    -----
    All hashes are content based. Covariance identity covers the complete native fit and relevant source/configuration, so touching timestamps cannot make a stale covariance valid. JSON publication uses a sibling temporary file and replacement; this is not a distributed-filesystem transaction.
    """
    root = Path(root)
    return dict(sources=source_hashes(root, repo), fitted_parameters=file_hash(root/'build/quadra_fitted_parameters.csv'))


def begin_covariance(root):
    """Save a pre-computation identity; invalidate any previous completed covariance stamp.

    Parameters
    ----------
    root
        Example directory used as the reference point for input/output paths.

    Notes
    -----
    All hashes are content based. Covariance identity covers the complete native fit and relevant source/configuration, so touching timestamps cannot make a stale covariance valid. JSON publication uses a sibling temporary file and replacement; this is not a distributed-filesystem transaction.
    """
    root = Path(root)
    # Remove the previous completion marker before starting new curvature. A
    # leftover covariance CSV alone can never certify a successful fresh run.
    (root/'build/covariance_provenance.json').unlink(missing_ok=True)
    atomic_json(root/'build/covariance_pending.json', covariance_inputs(root))


def finish_covariance(root):
    """Publish a covariance stamp only if inputs stayed unchanged during computation.

    Parameters
    ----------
    root
        Example directory used as the reference point for input/output paths.

    Notes
    -----
    All hashes are content based. Covariance identity covers the complete native fit and relevant source/configuration, so touching timestamps cannot make a stale covariance valid. JSON publication uses a sibling temporary file and replacement; this is not a distributed-filesystem transaction.
    """
    root = Path(root)
    pending = json.loads((root/'build/covariance_pending.json').read_text())
    if pending != covariance_inputs(root):
        raise ValueError('Inputs/source changed during curvature computation; rerun make uncertainty')
    names = ['coordinate_covariance.csv', 'uncertainty_coordinates.csv', 'hessian_coarse.csv', 'hessian_fine.csv', 'curvature_summary.json']
    atomic_json(root/'build/covariance_provenance.json', dict(schema=1, inputs=pending, artifacts={n:file_hash(root/'build'/n) for n in names}))
    (root/'build/covariance_pending.json').unlink()


def verify_covariance(root):
    """Fail closed on missing, stale, altered, or incomplete covariance provenance.

    Parameters
    ----------
    root
        Example directory used as the reference point for input/output paths.

    Returns
    -------
    Verified provenance dictionary; raises on an absent, stale, altered or incomplete stamp.

    Notes
    -----
    All hashes are content based. Covariance identity covers the complete native fit and relevant source/configuration, so touching timestamps cannot make a stale covariance valid. JSON publication uses a sibling temporary file and replacement; this is not a distributed-filesystem transaction.
    """
    root = Path(root)
    stamp = json.loads((root/'build/covariance_provenance.json').read_text())
    if stamp.get('schema') != 1 or stamp['inputs'] != covariance_inputs(root):
        raise ValueError('Covariance does not match source, configuration or fitted values; run make uncertainty')
    if set(stamp['artifacts']) != {'coordinate_covariance.csv','uncertainty_coordinates.csv','hessian_coarse.csv','hessian_fine.csv','curvature_summary.json'}:
        raise ValueError('Incomplete covariance provenance')
    for name, expected in stamp['artifacts'].items():
        if file_hash(root/'build'/name) != expected:
            raise ValueError('Covariance artifact changed: '+name)
    return stamp


def main():
    """Dispatch covariance begin/finish/verify commands from the example directory.

    Notes
    -----
    All hashes are content based. Covariance identity covers the complete native fit and relevant source/configuration, so touching timestamps cannot make a stale covariance valid. JSON publication uses a sibling temporary file and replacement; this is not a distributed-filesystem transaction.
    """
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['begin','finish','verify'])
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    {'begin':begin_covariance,'finish':finish_covariance,'verify':verify_covariance}[args.action](root)
    print('Covariance provenance: '+args.action+' passed')


if __name__ == '__main__':
    main()
