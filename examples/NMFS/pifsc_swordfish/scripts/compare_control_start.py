"""Check literal control seeds and compare independently fitted cold-start results."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from swordfish_inputs import load_inputs
from runtime import atomic_json
ROOT=Path(__file__).resolve().parents[1]


def parameters(path):
    """Read one complete native parameter CSV into a name/value map, rejecting duplicate names."""
    with path.open() as stream:
        rows=list(csv.DictReader(stream))
    result={r['parameter']:float(r['value']) for r in rows}
    if len(result)!=len(rows) or not all(math.isfinite(v) for v in result.values()):
        raise ValueError('Invalid parameter map')
    return result


def compare(strict,matched,output):
    """Validate both control-only starts; compare matched-fixed fits on an identical objective."""
    control=load_inputs(ROOT/'data/raw')['control']
    expected={r['name']:r['columns'][2] for r in control['parameters']}
    expected.update({f'RecruitmentDev_{y}':0. for y in range(1965,2022)})
    summaries={}
    for name,folder in [('strict',strict),('matched_fixed',matched)]:
        seed=parameters(folder/'initial_parameters.csv');target=expected.copy()
        if name=='matched_fixed':target['L_at_Amin_Mal_GP_1']=83.4
        if seed.keys()!=target.keys() or any(abs(seed[k]-v)>1e-12 for k,v in target.items()):
            raise ValueError('Starting values do not match the declared control experiment')
        summary=json.loads((folder/'summary.json').read_text())
        if not summary['converged'] or summary['gradient_norm']>1e-4 or summary['reads_SS3_estimates'] or summary['uses_SS3_covariance']:
            raise ValueError('Cold-start convergence or information-source check failed')
        summaries[name]=summary
    with (ROOT/'reference/fit_diagnostics/quadra_fit_summary.json').open() as stream:warm=json.load(stream)
    warm_parameters=parameters(ROOT/'reference/projection_uncertainty/quadra_fitted_parameters.csv')
    fitted=parameters(matched/'fitted_parameters.csv')
    errors=[]
    with (ROOT/'reference/fit_diagnostics/parameter_uncertainty.csv').open() as stream:
        for row in csv.DictReader(stream):
            name=row['parameter'];se=float(row['standard_error'])
            errors.append(dict(parameter=name,control_start=fitted[name],ss3_initialized=warm_parameters[name],difference=fitted[name]-warm_parameters[name],difference_in_reference_SE=(fitted[name]-warm_parameters[name])/se))
    maximum_standardized=max(abs(r['difference_in_reference_SE']) for r in errors)
    objective_difference=summaries['matched_fixed']['objective']-warm['objective']
    with (ROOT/'reference/projections/dynamic_unfished.csv').open() as stream:
        # The dynamic file's fished SSB is a saved canonical independent-fit output.
        ssb={int(r['year']):float(r['fished_spawning_biomass_tonnes']) for r in csv.DictReader(stream)}
    annual_errors=[]
    for label,folder in [('strict',strict),('matched_fixed',matched)]:
        with (folder/'quadra_fitted_annual.csv').open() as stream:
            for row in csv.DictReader(stream):
                year=int(row['year']);value=float(row['spawning_biomass'])
                annual_errors.append(dict(configuration=label,year=year,control_start_ssb_tonnes=value,ss3_initialized_ssb_tonnes=ssb[year],difference_tonnes=value-ssb[year]))
    matched_max=max(abs(r['difference_tonnes']) for r in annual_errors if r['configuration']=='matched_fixed')
    output.mkdir(parents=True,exist_ok=True)
    for filename,rows in [('parameter_comparison.csv',errors),('annual_ssb_comparison.csv',annual_errors)]:
        with (output/filename).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    result=dict(passed=abs(objective_difference)<1e-5 and maximum_standardized<.01 and matched_max<.1,
                strict=summaries['strict'],matched_fixed=summaries['matched_fixed'],ss3_initialized=warm,
                matched_objective_difference=objective_difference,maximum_parameter_difference_in_reference_SE=maximum_standardized,
                maximum_matched_annual_SSB_difference_tonnes=matched_max,
                interpretation='Strict controls use male L1=83.2, so they define a different objective from Report.sso. Matched-fixed uses 83.4 with every free coordinate still initialized from controls.',
                literal_start='The original full likelihood is not evaluable (empty predicted composition); the explicitly reported catch-only feasibility stage starts at the literal controls.',
                control_inputs_sha256={name:hashlib.sha256((ROOT/'data/raw'/name).read_bytes()).hexdigest() for name in ['swo2023_v004.dat','swo2023_v007.ctl']})
    atomic_json(output/'comparison_summary.json',result)
    if not result['passed']:raise ValueError('Matched-fixed cold fit has not recovered the reference optimum')
    print('Control starts verified; matched-fixed cold fit recovers the SS3-initialized optimum.')
    return result


def main():
    """Compare two completed independent cold runs without altering their optimization inputs."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strict',type=Path,required=True);parser.add_argument('--matched',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();compare(args.strict,args.matched,args.output)


if __name__=='__main__':main()
