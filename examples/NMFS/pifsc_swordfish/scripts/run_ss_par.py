"""Build and fit Quadra from a supplied SS3 .par vector in an isolated directory."""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from prepare_model import generate_header
from runtime import atomic_json, file_hash
from ss_par import parse
from swordfish_inputs import load_inputs

ROOT = Path(__file__).resolve().parents[1]


def run(parameter_file, output):
    """Validate/map .par values, compile against original dat/ctl, and fit from them.

    All generated fixtures and fitted outputs are isolated beneath a new output
    directory; the existing canonical fit is never replaced. The ordinary
    assessment executable reads .par-derived fixtures rather than Report.sso.
    Only the post-fit comparison reads the saved canonical Quadra summary.
    Compiler selection respects CXX; requires C++17 and the repository headers.
    """
    inputs = load_inputs(ROOT/'data/raw')
    parameter_hash = file_hash(parameter_file)
    native, metadata = parse(parameter_file, inputs)
    if file_hash(parameter_file) != parameter_hash:
        raise ValueError('Parameter file changed during parsing')
    output.mkdir(parents=True, exist_ok=False)
    scratch = output/'build'
    scratch.mkdir()
    identity = dict(parameter_file_sha256=parameter_hash,
                    observation_control_sha256=inputs['source_sha256'],
                    saved_ss3_metadata=metadata,
                    initialization='supplied ss.par, including fixed settings and historical deviations',
                    report_sso_used=False, saved_covariance_used=False)
    atomic_json(output/'input_manifest.json', identity)
    with (output/'initial_parameters.csv').open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['parameter', 'value'])
        writer.writerows(native.items())
    specs = inputs['control']['parameters']
    with (scratch/'biology_fitted_parameters.txt').open('w') as stream:
        stream.write(f'SWORDFISH_FITTED_PARAMETERS_V1\n{len(specs)}\n')
        for spec in specs:
            stream.write(f"{spec['name']} {native[spec['name']]:.17g}\n")
    with (scratch/'population_fitted_deviations.txt').open('w') as stream:
        stream.write('SWORDFISH_DEVIATIONS_V1\n57\n')
        for year in range(1965, 2022):
            stream.write(f"{year} {native[f'RecruitmentDev_{year}']:.17g}\n")
    generate_header(inputs, scratch/'swordfish_input.hpp')
    repository = ROOT.parents[2]
    binary = output/'swordfish_ss_par'
    command = shlex.split(os.environ.get('CXX', 'c++')) + [
        '-std=c++17', '-O2', '-I'+str(scratch), '-I'+str(ROOT/'include'),
        '-I'+str(repository), '-I'+str(repository/'external/eigen'),
        str(ROOT/'src/swordfish_assessment.cpp'), '-o', str(binary)]
    print('Run directory: '+str(output), flush=True)
    print('Building Quadra from original dat/ctl and supplied .par values...', flush=True)
    subprocess.run(command, check=True)
    gradient = subprocess.run([str(binary), '--gradient'], cwd=output,
                              text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output/'initial_gradient.log').write_text(gradient.stdout)
    print(gradient.stdout, end='', flush=True)
    gradient.check_returncode()
    match = re.search(r'Objective: (\S+); gradient norm: (\S+);', gradient.stdout)
    if not match:
        raise ValueError('Missing initial gradient result')
    initial_objective, initial_gradient = map(float, match.groups())
    with (output/'fit.log').open('w') as log:
        with subprocess.Popen([str(binary), '--fit', '200'], cwd=output, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1) as process:
            for line in process.stdout:
                sys.stdout.write(line); sys.stdout.flush(); log.write(line); log.flush()
            status = process.wait()
            if status:
                raise subprocess.CalledProcessError(status, process.args)
    if file_hash(parameter_file) != identity['parameter_file_sha256'] or load_inputs(ROOT/'data/raw')['source_sha256'] != inputs['source_sha256']:
        raise ValueError('Input files changed during the fit')
    fitted = json.loads((scratch/'quadra_fit_summary.json').read_text())
    reference = json.loads((ROOT/'reference/fit_diagnostics/quadra_fit_summary.json').read_text())
    ss3_difference = fitted['objective']-metadata['ss3_objective']
    reference_difference = fitted['objective']-reference['objective']
    with (scratch/'quadra_fitted_parameters.csv').open() as stream:
        fit_values = {r['parameter']: float(r['value']) for r in csv.DictReader(stream)}
    with (ROOT/'reference/projection_uncertainty/quadra_fitted_parameters.csv').open() as stream:
        reference_values = {r['parameter']: float(r['value']) for r in csv.DictReader(stream)}
    with (ROOT/'reference/fit_diagnostics/parameter_uncertainty.csv').open() as stream:
        parameter_rows = [dict(parameter=r['parameter'], ss_par_initialized=fit_values[r['parameter']],
                               existing_quadra=reference_values[r['parameter']],
                               difference=fit_values[r['parameter']]-reference_values[r['parameter']],
                               difference_in_reference_SE=(fit_values[r['parameter']]-reference_values[r['parameter']])/float(r['standard_error']))
                          for r in csv.DictReader(stream)]
    with (ROOT/'reference/projections/dynamic_unfished.csv').open() as stream:
        reference_ssb = {int(r['year']): float(r['fished_spawning_biomass_tonnes']) for r in csv.DictReader(stream)}
    with (scratch/'quadra_fitted_annual.csv').open() as stream:
        annual_rows = [dict(year=int(r['year']), ss_par_initialized_ssb_tonnes=float(r['spawning_biomass']),
                            existing_quadra_ssb_tonnes=reference_ssb[int(r['year'])],
                            difference_tonnes=float(r['spawning_biomass'])-reference_ssb[int(r['year'])]) for r in csv.DictReader(stream)]
    for filename, rows in [('parameter_comparison.csv', parameter_rows), ('annual_ssb_comparison.csv', annual_rows)]:
        with (output/filename).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
    result = dict(initial_objective=initial_objective, initial_gradient_norm=initial_gradient,
                  quadra=fitted, ss3_objective=metadata['ss3_objective'],
                  quadra_minus_ss3_objective=ss3_difference,
                  quadra_minus_existing_quadra_objective=reference_difference,
                  objective_tolerance=1e-4,
                  passed=fitted['converged'] and fitted['gradient_norm'] <= 1e-4
                  and abs(ss3_difference) < 1e-4 and abs(reference_difference) < 1e-4,
                  interpretation='Agreement from a saved SS3 fitted vector; no SS3 rerun or original pre-fit starting vector is established.',
                  profiled_q='Supplied Q values are mapped; Quadra analytically profiles Q at every objective evaluation.',
                  main_deviation_constraint='Final main deviation is reconstructed as negative sum; input sum tolerance 1e-8.',
                  binary_sha256=file_hash(binary))
    result['maximum_parameter_difference_in_reference_SE'] = max(abs(r['difference_in_reference_SE']) for r in parameter_rows)
    result['maximum_annual_ssb_difference_tonnes'] = max(abs(r['difference_tonnes']) for r in annual_rows)
    atomic_json(output/'comparison_summary.json', result)
    if not result['passed']:
        raise ValueError('SS3 .par agreement check failed; inspect comparison_summary.json')
    print('PASS: Quadra fit agrees with the SS3 .par objective and existing Quadra optimum.', flush=True)
    return result


def main():
    """Accept a local .par file and optional new output directory; show optimizer progress."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parameter_file', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = args.output or ROOT/'build/ss_par'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    run(args.parameter_file.resolve(), output.resolve())


if __name__ == '__main__':
    main()
