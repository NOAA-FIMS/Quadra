"""Clean native-only swordfish release: frozen sources, local starts, uncertainty and PDF."""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import random
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from native_input import load, encode
from runtime import atomic_json, file_hash
from compare_benchmarks import compare as compare_points

ROOT=Path(__file__).resolve().parents[1]
TARGETS=['ss_par/initial_parameters.csv','ss_par/input_manifest.json',
         'fit_diagnostics/quadra_fit_summary.json','fit_diagnostics/benchmark_comparison.json',
         'fit_diagnostics/uncertainty_comparison.json']


def hashes(directory):
    """Hash a source/artifact tree in relative-path order, excluding Python caches."""
    return {str(p.relative_to(directory)):file_hash(p) for p in sorted(directory.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts and p.name!='.DS_Store'}


def copy_snapshot(case, bundle):
    """Copy model/framework source and selected comparison targets, never raw SS3 files or binaries."""
    repo=case.parents[2];workspace=bundle/'workspace'
    ignore=shutil.ignore_patterns('build','__pycache__','.DS_Store','.git')
    for name in ['core','model','math','include','external/eigen','external/LBFGSpp']:
        shutil.copytree(repo/name,workspace/name,ignore=ignore)
    (workspace/'tests').mkdir()
    for p in (repo/'tests').glob('*.cpp'):shutil.copyfile(p,workspace/'tests'/p.name)
    shutil.copyfile(repo/'Makefile',workspace/'Makefile')
    frozen=workspace/'examples/NMFS/pifsc_swordfish';frozen.mkdir(parents=True)
    for name in ['include','src','scripts','tests']:shutil.copytree(case/name,frozen/name,ignore=ignore)
    for name in ['Makefile','swordfish']:shutil.copy2(case/name,frozen/name)
    for name in TARGETS:
        target=frozen/'reference'/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(case/'reference'/name,target)
    return frozen


def run_stage(name, command, cwd, bundle, records, timeout=1800):
    """Stream stage output to console/log, enforce a deadline and terminate its process group on failure."""
    log=bundle/'logs'/f'{len(records)+1:02d}_{name}.log';log.parent.mkdir(exist_ok=True)
    print('Running '+name+'; log '+str(log),flush=True)
    started=time.monotonic();errors=[]
    with log.open('w') as stream:
        process=subprocess.Popen(command,cwd=cwd,text=True,stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,start_new_session=True)
        def pump():
            """Drain child output continuously, retaining logs if console output closes."""
            console=True
            try:
                for line in process.stdout:
                    stream.write(line);stream.flush()
                    if console:
                        try:sys.stdout.write(line);sys.stdout.flush()
                        except (BrokenPipeError,OSError):console=False
            except BaseException as error:
                errors.append(error)
                try:os.killpg(process.pid,signal.SIGTERM)
                except ProcessLookupError:pass
        reader=threading.Thread(target=pump,daemon=True);reader.start()
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
            raise
        finally:
            reader.join(timeout=10)
            if not reader.is_alive():process.stdout.close()
            records.append(dict(name=name,command=[str(v) for v in command],seconds=time.monotonic()-started,
                                returncode=process.returncode,log=str(log.relative_to(bundle))))
        if reader.is_alive() or errors:raise RuntimeError('Stage log streaming failed: '+name)
    if process.returncode:raise RuntimeError(f'{name} failed ({process.returncode}); see {log.name}')


def parameter_values(path):
    """Read a finite, duplicate-free complete native parameter map for seed comparison."""
    with path.open() as stream:rows=list(csv.DictReader(stream))
    values={r['parameter']:float(r['value']) for r in rows}
    if len(values)!=len(rows) or not all(math.isfinite(v) for v in values.values()):raise ValueError('Invalid seed parameters')
    return values


def verify_seed(document, expected):
    """Require the release baseline to use the preserved .par-imported vector before Quadra fitting."""
    values={p['name']:p['value'] for p in document['parameters']}
    values.update({f"RecruitmentDev_{d['year']}":d['value'] for d in document['recruitment_deviations']})
    if values.keys()!=expected.keys() or any(abs(values[k]-expected[k])>1e-12 for k in expected):
        raise ValueError('Release baseline is not the recorded ss.par initial vector')


def coordinate_schema(document):
    """Match the model's independent bounded-logit coordinates, excluding the final main deviation."""
    schema=[]
    for p in document['parameters']:
        if p['phase']>0:schema.append((p['name'],p['lower'],p['upper'],p['value']))
    rec=document['model']['recruitment'];deviations={r['year']:r['value'] for r in document['recruitment_deviations']}
    for year in range(rec['early_start'],rec['main_end']):
        schema.append((f'RecruitmentDev_{year}',rec['minimum_deviation'],rec['maximum_deviation'],deviations[year]))
    return schema


def metric_factor(path, names):
    """Read the fresh labeled covariance and compute an unregularized numeric Cholesky factor."""
    with path.open() as stream:rows=list(csv.reader(stream))
    if rows[0]!=['quantity']+names or len(rows)!=len(names)+1:raise ValueError('Perturbation covariance labels/dimensions')
    matrix=[]
    for name,row in zip(names,rows[1:]):
        if row[0]!=name or len(row)!=len(names)+1:raise ValueError('Perturbation covariance row order')
        matrix.append([float(v) for v in row[1:]])
    if not all(math.isfinite(v) for row in matrix for v in row):raise ValueError('Nonfinite perturbation covariance')
    n=len(names);factor=[[0.]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            if abs(matrix[i][j]-matrix[j][i])>1e-10*max(1.,abs(matrix[i][j])):raise ValueError('Asymmetric perturbation covariance')
            value=matrix[i][j]-sum(factor[i][k]*factor[j][k] for k in range(j))
            if i==j:
                if value<=0:raise ValueError('Perturbation covariance is not positive definite')
                factor[i][j]=math.sqrt(value)
            else:factor[i][j]=value/factor[j][j]
    return factor


def perturbed_input(document, factor, seed, scale):
    """Perturb the .par starting coordinates along fresh covariance directions; preserve fixed settings."""
    result=copy.deepcopy(document);schema=coordinate_schema(document);rng=random.Random(seed)
    z=[rng.gauss(0,scale) for _ in schema];values={}
    for i,(name,lo,hi,value) in enumerate(schema):
        x=math.log((value-lo)/(hi-value))+sum(factor[i][j]*z[j] for j in range(len(z)))
        probability=1/(1+math.exp(-x)) if x>=0 else math.exp(x)/(1+math.exp(x))
        values[name]=lo+(hi-lo)*probability
    for p in result['parameters']:
        if p['phase']>0:p['value']=values[p['name']]
    rec=result['model']['recruitment']
    for d in result['recruitment_deviations']:
        if d['year']<rec['main_end']:d['value']=values[f"RecruitmentDev_{d['year']}"]
    final=-sum(d['value'] for d in result['recruitment_deviations'] if rec['main_start']<=d['year']<rec['main_end'])
    next(d for d in result['recruitment_deviations'] if d['year']==rec['main_end'])['value']=final
    encode(result);return result


def accept(bundle, frozen, records):
    """Gate convergence, SS3 points/SEs, start recovery and run-bound input/binary/covariance hashes."""
    baseline=bundle/'results/baseline';fit=json.loads((baseline/'quadra_fit_summary.json').read_text())
    uncertainty=json.loads((baseline/'reference_uncertainty_summary.json').read_text())
    if not fit['converged'] or fit['gradient_norm']>1e-4 or not uncertainty['passed']:raise ValueError('Baseline convergence/uncertainty gate')
    expected=json.loads((frozen/'reference/ss_par/input_manifest.json').read_text())['saved_ss3_metadata']['ss3_objective']
    canonical=json.loads((frozen/'reference/fit_diagnostics/quadra_fit_summary.json').read_text())['objective']
    if abs(fit['objective']-expected)>1e-4 or abs(fit['objective']-canonical)>1e-5:raise ValueError('Baseline objective comparison failed')
    points=json.loads((frozen/'reference/fit_diagnostics/benchmark_comparison.json').read_text())
    point_check=compare_points(baseline/'management_reference_points.csv',{r['quantity']:r['ss3'] for r in points['checks']})
    if not point_check['passed']:raise ValueError('SS3 reference-point gate')
    atomic_json(bundle/'reference_point_comparison.json',point_check)
    expected_se={r['quantity']:r['ss3_se'] for r in json.loads((frozen/'reference/fit_diagnostics/uncertainty_comparison.json').read_text())['comparisons']}
    prefix={'SSB_tonnes':'SSB','dead_catch_tonnes':'Dead_Catch','annual_F_SS3_units4_ages1_10':'annF'};suffix={'unfished':'unfished','Btgt20':'Btgt','SPR20':'SPR','MSY':'MSY'};checks=[]
    with (baseline/'management_reference_point_uncertainty.csv').open() as stream:
        for row in csv.DictReader(stream):
            if row['quantity'] in prefix:
                name=prefix[row['quantity']]+'_'+suffix[row['benchmark']]
                if name in expected_se:
                    actual=float(row['standard_error']);expected_value=expected_se[name]
                    checks.append(dict(quantity=name,quadra_SE=actual,ss3_SE=expected_value,relative_difference=abs(actual/expected_value-1)))
    if len(checks)!=10 or any(r['relative_difference']>=.01 for r in checks):raise ValueError('SS3 reference SE gate')
    atomic_json(bundle/'reference_SE_comparison.json',dict(passed=True,tolerance_relative=.01,checks=checks))
    with (baseline/'objective_gradient.csv').open() as stream:center=list(csv.DictReader(stream))
    names=[r['parameter'] for r in center];factor=metric_factor(baseline/'reference_coordinate_covariance.csv',names)
    deviations=[math.sqrt(sum(v*v for v in row)) for row in factor]
    recoveries=[];binary_hash=file_hash(frozen/'build/swordfish_native')
    for info in records:
        folder=bundle/'results'/info['name'];status=json.loads((folder/'run_status.json').read_text())
        if status['status']!='completed' or status['binary_sha256']!=binary_hash or status['input_sha256']!=file_hash(bundle/'inputs'/info['input']):raise ValueError('Run/input/binary identity failed')
        if info['name']=='baseline':continue
        if status['optimizer_metric_sha256']!=file_hash(baseline/'reference_coordinate_covariance.csv') or not status['optimizer_metric_unchanged']:raise ValueError('Perturbation metric identity failed')
        result=json.loads((folder/'quadra_fit_summary.json').read_text())
        with (folder/'objective_gradient.csv').open() as stream:coordinates=list(csv.DictReader(stream))
        if [r['parameter'] for r in coordinates]!=names:raise ValueError('Recovery coordinate order')
        distance=max(abs(float(a['unconstrained_value'])-float(b['unconstrained_value']))/sd for a,b,sd in zip(coordinates,center,deviations))
        if not result['converged'] or result['gradient_norm']>1e-4 or abs(result['objective']-fit['objective'])>1e-5 or distance>.01:raise ValueError('Perturbed start did not recover the optimum')
        recoveries.append(dict(**info,initial_objective=result['initial_objective'],objective=result['objective'],gradient_norm=result['gradient_norm'],maximum_coordinate_difference_in_SE=distance,optimizer_iterations=result['optimizer_iterations']))
    return dict(passed=True,baseline=fit,ss3_objective=expected,ss3_objective_difference=fit['objective']-expected,
                reference_points_passed=True,reference_SEs_passed=True,uncertainty=uncertainty,perturbed_starts=recoveries,
                interpretation='Native-input reproducibility and local start recovery; not global optimality or bootstrap/projection validation.')


def release(case, input_path, destination, pdf_python):
    """Freeze, build, fit and publish a new local native bundle only after every gate and PDF succeeds."""
    case=Path(case).resolve();destination=Path(destination).absolute()
    if destination.exists() or destination.is_symlink():raise FileExistsError('Use a new release destination')
    source=input_path.read_bytes();document,_=load(source.decode('utf-8'))
    verify_seed(document,parameter_values(case/'reference/ss_par/initial_parameters.csv'))
    compiler=os.environ.get('CXX','c++')
    if not shutil.which(compiler):raise RuntimeError('A C++17 compiler is required')
    try:reportlab_version=subprocess.check_output([str(pdf_python),'-c','import reportlab; print(reportlab.Version)'],text=True).strip()
    except (subprocess.CalledProcessError,OSError) as error:raise RuntimeError('PDF generation requires ReportLab; select its Python with --pdf-python') from error
    destination.parent.mkdir(parents=True,exist_ok=True);lock=destination.with_name('.'+destination.name+'.lock');lock.mkdir()
    staging=Path(tempfile.mkdtemp(prefix='.'+destination.name+'.inprogress-',dir=destination.parent));stages=[]
    try:
        frozen=copy_snapshot(case,staging);sources=hashes(staging/'workspace')
        if any(file_hash(case.parents[2]/name)!=digest for name,digest in sources.items()):raise ValueError('Source changed while copying snapshot')
        atomic_json(staging/'source_snapshot.json',dict(files=sources,raw_ss3_files_included=False))
        if (frozen/'data/raw').exists() or list((frozen/'build').glob('*')):raise ValueError('Snapshot is not a clean native build')
        (staging/'inputs').mkdir();(staging/'inputs/model.json').write_bytes(source)
        (staging/'results').mkdir();make=['make','CXX='+compiler,'PYTHON='+sys.executable]
        run_stage('clean_native_build',make+['native-swordfish','test-portable','build/test_optimizer_scaling'],frozen,staging,stages)
        run_stage('optimizer_scaling_checks',['./build/test_optimizer_scaling'],frozen,staging,stages)
        run_stage('native_schema_checks',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_native_input.py'],frozen,staging,stages)
        run_stage('release_guardrails',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_native_release.py'],frozen,staging,stages)
        run_stage('framework_regressions',make+['test-framework'],frozen,staging,stages)
        run_stage('baseline_fit_uncertainty',[sys.executable,'swordfish','--input',str(staging/'inputs/model.json'),'--output',str(staging/'results/baseline'),'--uncertainty'],frozen,staging,stages)
        records=[dict(name='baseline',input='model.json')]
        covariance=staging/'results/baseline/reference_coordinate_covariance.csv'
        schema=coordinate_schema(document);factor=metric_factor(covariance,[r[0] for r in schema])
        for index,(seed,scale) in enumerate([(101,.1),(202,.25),(303,.5)],1):
            name=f'start-{index:02d}';filename=name+'.json';candidate=perturbed_input(document,factor,seed,scale)
            atomic_json(staging/'inputs'/filename,candidate)
            run_stage(name,[sys.executable,'swordfish','--input',str(staging/'inputs'/filename),'--output',str(staging/'results'/name),'--max-iterations','300','--optimizer-metric',str(covariance)],frozen,staging,stages)
            records.append(dict(name=name,input=filename,seed=seed,scale=scale))
        acceptance=accept(staging,frozen,records);atomic_json(staging/'acceptance.json',acceptance);atomic_json(staging/'stages.json',stages)
        run_stage('summary_pdf',[str(pdf_python),'scripts/build_native_release_report.py',str(staging)],frozen,staging,stages)
        current={k:v for k,v in hashes(staging/'workspace').items() if 'build' not in Path(k).parts}
        if current!=sources:raise ValueError('Frozen source or comparison targets changed')
        atomic_json(staging/'stages.json',stages)
        (staging/'README.md').write_text('# Native swordfish release bundle\n\nSee acceptance.json, release_manifest.json, output/pdf/native_release_summary.pdf and results/baseline/.\nSources are frozen in workspace/; input JSON and all seeded perturbations are under inputs/.\nThe baseline covariance scales optimizer steps for local-start checks only.\nThis local bundle contains observations and should be shared deliberately.\n')
        manifest=dict(schema=1,status='passed',created_utc=datetime.now(timezone.utc).isoformat(),sources=sources,
                      input_sha256=file_hash(staging/'inputs/model.json'),stages=stages,acceptance=acceptance,
                      environment=dict(python=sys.version,pdf_python=str(pdf_python),reportlab=reportlab_version,platform=platform.platform(),compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]),
                      scope='Clean native fit, local perturbed starts and reference-point uncertainty; no raw SS3 files used.')
        manifest['artifacts']=hashes(staging);atomic_json(staging/'release_manifest.json',manifest)
        if destination.exists() or destination.is_symlink():raise FileExistsError('Destination appeared; refusing replacement')
        os.rename(staging,destination);print('Published native release: '+str(destination),flush=True);return destination
    except BaseException as error:
        if staging.exists():
            atomic_json(staging/'failure.json',dict(status='failed',error=str(error),stages=stages))
            failed=staging.with_name(staging.name.replace('.inprogress-','.failed-'));os.rename(staging,failed)
            print('Failed native run retained at '+str(failed),file=sys.stderr)
        raise
    finally:lock.rmdir()


def main():
    """Run a clean native-only release using the recorded .par-derived starting JSON."""
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path);parser.add_argument('--pdf-python',type=Path,default=Path(sys.executable))
    args=parser.parse_args();output=args.output or ROOT/'runs'/('native-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    release(ROOT,args.input.resolve(),output,args.pdf_python.absolute())


if __name__=='__main__':main()
