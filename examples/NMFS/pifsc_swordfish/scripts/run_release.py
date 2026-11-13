"""Run a clean swordfish release acceptance pipeline and atomically publish a new bundle."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from runtime import atomic_json, file_hash, source_hashes


def verify_sources(case):
    """Require all eight original assessment files to match the recorded source manifest.

    Parameters
    ----------
    case
        Example directory containing source, original data and build outputs.

    Notes
    -----
    The supported unattended workflow uses a frozen source/input copy and serial stages. Scratch outputs remain unpublished until every numerical gate passes. A destination lock coordinates this runner only; an uncatchable termination can leave a lock requiring operator inspection. Original data and failed runs stay local.
    """
    entries = json.loads((case/'reference/source_manifest.json').read_text())
    if len(entries) != 8:
        raise ValueError('Expected eight source manifest entries')
    for row in entries:
        if Path(row['file']).name != row['file'] or file_hash(case/'data/raw'/row['file']) != row['sha256']:
            raise ValueError('Source manifest mismatch: '+row['file'])


def copy_snapshot(case, destination):
    """Copy case and bundled framework source, excluding existing binaries and run caches.

    Parameters
    ----------
    case
        Example directory containing source, original data and build outputs.
    destination
        New output/snapshot directory; completed release destinations are never overwritten.

    Returns
    -------
    Path to the isolated example inside the copied workspace.

    Notes
    -----
    The supported unattended workflow uses a frozen source/input copy and serial stages. Scratch outputs remain unpublished until every numerical gate passes. A destination lock coordinates this runner only; an uncatchable termination can leave a lock requiring operator inspection. Original data and failed runs stay local.
    """
    repo = case.parents[2]
    frozen = destination/'workspace'
    for name in ['core','model','math','include','external/eigen','external/LBFGSpp']:
        shutil.copytree(repo/name, frozen/name)
    (frozen/'tests').mkdir()
    for path in (repo/'tests').glob('*.cpp'):
        shutil.copyfile(path, frozen/'tests'/path.name)
    shutil.copyfile(repo/'Makefile', frozen/'Makefile')
    target = frozen/'examples/NMFS/pifsc_swordfish'
    shutil.copytree(case, target, ignore=shutil.ignore_patterns('build','runs','__pycache__','.DS_Store','output','tmp'))
    return target


def run_stage(name, command, cwd, bundle, records, timeout=1800):
    """Execute one checked stage with a 30-minute limit; terminate its process group on interruption.

    Parameters
    ----------
    name
        Exact section, output, field or stage identifier; matching is case sensitive.
    command
        Argument list passed directly to a subprocess, never interpolated into a shell command.
    cwd
        Working directory for the child stage.
    bundle
        In-progress bundle directory receiving stage logs.
    records
        Mutable list to which execution status, log path and elapsed time are appended.
    timeout
        Maximum stage runtime in seconds; timeout triggers process-group termination.

    Notes
    -----
    The supported unattended workflow uses a frozen source/input copy and serial stages. Scratch outputs remain unpublished until every numerical gate passes. A destination lock coordinates this runner only; an uncatchable termination can leave a lock requiring operator inspection. Original data and failed runs stay local.
    """
    log = bundle/'logs'/f'{len(records)+1:02d}_{name}.log'
    log.parent.mkdir(exist_ok=True)
    print('Running '+name, flush=True)
    start = time.monotonic()
    with log.open('w') as stream:
        process = subprocess.Popen(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            process.wait(timeout=timeout)
        except BaseException:
            try:os.killpg(process.pid,signal.SIGTERM)
            except ProcessLookupError:pass
            try:process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                process.wait()
            records.append(dict(name=name,command=command,seconds=time.monotonic()-start,returncode=process.returncode,log=str(log.relative_to(bundle)),interrupted=True))
            raise
    records.append(dict(name=name, command=command, seconds=time.monotonic()-start, returncode=process.returncode, log=str(log.relative_to(bundle))))
    if process.returncode:
        raise RuntimeError(f'{name} failed ({process.returncode}); see {log.name}')


def accept_results(case):
    """Check machine-readable scientific gates, fit convergence and covariance identity.

    Parameters
    ----------
    case
        Example directory containing source, original data and build outputs.

    Returns
    -------
    Parsed machine-readable numerical gates, all required to pass.

    Notes
    -----
    The supported unattended workflow uses a frozen source/input copy and serial stages. Scratch outputs remain unpublished until every numerical gate passes. A destination lock coordinates this runner only; an uncatchable termination can leave a lock requiring operator inspection. Original data and failed runs stay local.
    """
    from runtime import verify_covariance
    verify_covariance(case)
    gates = ['quadra_fit_summary.json','quadra_fit_comparison.json','benchmark_comparison.json','robustness_summary.json','uncertainty_comparison.json','projections/dynamic_unfished_comparison.json','projections/historical_SPR_comparison.json','projections/projection_uncertainty_validation.json']
    results = {}
    for name in gates:
        value = json.loads((case/'build'/name).read_text())
        passed=value.get('converged',value.get('passed'))
        if name=='quadra_fit_comparison.json':passed=all(row.get('passed') is True for row in value.values())
        if passed is not True:
            raise ValueError('Acceptance gate failed: '+name)
        results[name] = value
    return results


def release(case, destination):
    """Lock a new destination, test an immutable copy, and publish only after all gates pass.

    Parameters
    ----------
    case
        Example directory containing source, original data and build outputs.
    destination
        New output/snapshot directory; completed release destinations are never overwritten.

    Returns
    -------
    Published destination path after every stage and gate passes.

    Notes
    -----
    The supported unattended workflow uses a frozen source/input copy and serial stages. Scratch outputs remain unpublished until every numerical gate passes. A destination lock coordinates this runner only; an uncatchable termination can leave a lock requiring operator inspection. Original data and failed runs stay local.
    """
    case, destination = Path(case).resolve(), Path(destination).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError('Destination already exists; use a new run name')
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock = destination.with_name('.'+destination.name+'.lock')
    # Reserving the destination name does not create a public result directory.
    # Publication happens only after every gate; crashes preserve private staging.
    lock.mkdir()  # Concurrent runs targeting the same name fail before changing outputs.
    staging = Path(tempfile.mkdtemp(prefix='.'+destination.name+'.inprogress-',dir=destination.parent))
    records = []
    try:
        verify_sources(case)
        initial = source_hashes(case)
        frozen = copy_snapshot(case,staging)
        if source_hashes(frozen) != initial or source_hashes(case) != initial:
            raise ValueError('Source changed while the snapshot was copied')
        compiler = os.environ.get('CXX','c++')
        if not shutil.which(compiler) or not shutil.which('Rscript'):
            raise RuntimeError('Release requires a C++17 compiler and Rscript on PATH')
        environment = dict(python=sys.version, platform=platform.platform(), compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0], R=subprocess.check_output(['Rscript','--version'],stderr=subprocess.STDOUT,text=True).strip())
        make = ['make','CXX='+compiler,'PYTHON='+sys.executable]
        run_stage('clean_build_all',make+['build-all'],frozen,staging,records)
        run_stage('unit_and_ss3_audits',make+['test'],frozen,staging,records)
        run_stage('framework_regressions',make+['test-framework'],frozen,staging,records)
        run_stage('build_fit',make+['build/swordfish_assessment'],frozen,staging,records)
        run_stage('fit',[str(frozen/'build/swordfish_assessment'),'--fit','500'],frozen,staging,records)
        run_stage('fit_and_benchmarks',make+['compare-fit','benchmarks'],frozen,staging,records)
        run_stage('fresh_curvature',make+['uncertainty'],frozen,staging,records)
        run_stage('perturbed_starts',make+['robustness'],frozen,staging,records)
        run_stage('projection_acceptance',make+['test-projection','test-projection-uncertainty'],frozen,staging,records)
        run_stage('figures',make+['plot-projections','plot-projection-uncertainty'],frozen,staging,records)
        gates = accept_results(frozen)
        if source_hashes(frozen) != initial:
            raise ValueError('Frozen source/input changed during release')
        artifacts = {str(p.relative_to(staging)):file_hash(p) for p in (frozen/'build').rglob('*') if p.is_file()}
        atomic_json(staging/'release_manifest.json',dict(schema=1,status='passed',created_utc=datetime.now(timezone.utc).isoformat(),environment=environment,sources=initial,stages=records,gates=gates,artifacts=artifacts,scope='operational acceptance for the supplied swordfish case; independent scientific approval and published bootstrap replay remain outstanding'))
        if destination.exists() or destination.is_symlink():
            raise FileExistsError('Destination appeared during the run; refusing replacement')
        # Staging and destination share a parent/filesystem. This rename is the
        # publication boundary; direct scratch writers are not covered by it.
        os.rename(staging,destination)
        print('Published '+str(destination),flush=True)
        return destination
    except BaseException as error:
        if staging.exists():
            atomic_json(staging/'failure.json',dict(status='failed',error=str(error),stages=records))
            failed=staging.with_name(staging.name.replace('.inprogress-','.failed-'))
            os.rename(staging,failed)
            print('Failed run retained at '+str(failed),file=sys.stderr)
        raise
    finally:
        lock.rmdir()


def main():
    """Parse a non-overwriting output path and run release acceptance.

    Notes
    -----
    The supported unattended workflow uses a frozen source/input copy and serial stages. Scratch outputs remain unpublished until every numerical gate passes. A destination lock coordinates this runner only; an uncatchable termination can leave a lock requiring operator inspection. Original data and failed runs stay local.
    """
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='New destination; default runs/swordfish-UTC timestamp')
    args=parser.parse_args()
    case=Path(__file__).resolve().parents[1]
    output=args.output or case/'runs'/('swordfish-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    release(case,output)


if __name__=='__main__':
    main()
