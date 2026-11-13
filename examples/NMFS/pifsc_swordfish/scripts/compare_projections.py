"""Validate dynamic unfished reproduction and describe, without calibrating to, published projections."""
import csv
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'build/projections'

def read(name):
    """Read a named projection CSV from the example build output folder.

    Parameters
    ----------
    name
        Exact section, output, field or stage identifier; matching is case sensitive.

    Notes
    -----
    Dynamic unfished and historical SPR values have numerical pass/fail tolerances. Published projection medians are descriptive comparisons, not equivalent validation targets: the original bootstrap models and forecast settings were not supplied.
    """
    with (OUT/name).open() as f:
        return list(csv.DictReader(f))

def main():
    """Execute the compare projections command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    Dynamic unfished and historical SPR values have numerical pass/fail tolerances. Published projection medians are descriptive comparisons, not equivalent validation targets: the original bootstrap models and forecast settings were not supplied.
    """
    lines = (ROOT/'data/raw/Report.sso').read_text().splitlines()
    start = lines.index('Spawning_Biomass_Report_1 No_fishery_for_Z=M_and_dynamic_Bzero')
    ssb = {}
    for line in lines[start+3:]:
        cells = line.split()
        if not cells:
            break
        if len(cells) == 3 and cells[1] == 'TIME':
            ssb[int(cells[0])] = float(cells[2])
    start = lines.index('NUMBERS_AT_AGE_Annual_1 No_fishery_for_Z=M_and_dynamic_Bzero')
    numbers = {}
    for line in lines[start+2:]:
        cells = line.split()
        if not cells:
            break
        if len(cells) != 19:
            raise ValueError('Unrecognized SS3 unfished numbers table')
        if 1975 <= int(cells[2]) <= 2021:
            for age, value in enumerate(cells[3:]):
                numbers[int(cells[2]), int(cells[1]), age] = float(value)
    spr_start = next(i for i, line in enumerate(lines) if line.startswith('Yr Era Bio_all Bio_Smry SSBzero'))
    spr_targets = {}
    for line in lines[spr_start+1:]:
        cells = line.split()
        if not cells:
            break
        if len(cells) >= 8 and cells[0].isdigit() and cells[1] == 'TIME' and 1975 <= int(cells[0]) <= 2021:
            spr_targets[int(cells[0])] = float(cells[7])
    spr_checks = [dict(year=int(r['year']), quadra=float(r['raw_SPR']), ss3=spr_targets[int(r['year'])]) for r in read('historical_SPR.csv')]
    for row in spr_checks:
        row['relative_error'] = abs(row['quadra']/row['ss3']-1)
    spr_result = dict(passed=len(spr_checks) == 47 and all(r['relative_error'] < 2e-5 for r in spr_checks), comparisons=len(spr_checks), maximum_relative_error=max(r['relative_error'] for r in spr_checks), tolerance_relative=2e-5, note='Raw SPR uses each historical year selectivity and fishing pattern; report SPR_report is 1-SPR.')
    (OUT/'historical_SPR_comparison.json').write_text(json.dumps(spr_result, indent=2)+'\n')
    history = read('dynamic_unfished.csv')
    checks = [dict(quantity='SSB', year=int(r['year']), quadra=float(r['spawning_biomass_tonnes']), ss3=ssb[int(r['year'])]) for r in history]
    for row in read('dynamic_unfished_numbers.csv'):
        key = tuple(int(row[k]) for k in ['year', 'sex', 'age'])
        checks.append(dict(quantity=f'numbers_sex{key[1]}_age{key[2]}', year=key[0], quadra=float(row['numbers_SS3_annual_thousands']), ss3=numbers[key]))
    for check in checks:
        check['relative_error'] = abs(check['quadra']-check['ss3'])/max(1., abs(check['ss3']))
        check['passed'] = check['relative_error'] < 2e-5
    mean = sum(float(r['spawning_biomass_tonnes']) for r in history if int(r['year']) >= 2017)/5
    report_mean = sum(ssb[y] for y in range(2017, 2022))/5
    dynamic = dict(passed=len(checks) == 1551 and all(c['passed'] for c in checks), comparisons=len(checks), tolerance_relative_with_unit_floor=2e-5, maximum_relative_error=max(c['relative_error'] for c in checks), mean_2017_2021_tonnes=mean, rounded_ss3_report_mean_tonnes=report_mean, published_assessment_mean_tonnes=95732, difference_from_published_tonnes=mean-95732, twenty_percent_dynamic_target_tonnes=.2*mean, reason_for_separate_baselines='Supplied Report.sso is a rerun. Published PDF mean and rounded report mean differ; retain both.')
    (OUT/'dynamic_unfished_comparison.json').write_text(json.dumps(dynamic, indent=2)+'\n')
    with (OUT/'dynamic_unfished_comparison.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=checks[0]);writer.writeheader();writer.writerows(checks)
    projections = {(int(r['scenario']), int(r['year'])): r for r in read('deterministic_projections.csv')}
    assert len(projections) == 50
    with (ROOT/'reference/published_projections.csv').open() as f:
        published = list(csv.DictReader(f))
    comparisons = []
    for r in published:
        key = int(r['scenario']), int(r['year']);field = 'spawning_biomass_tonnes' if r['quantity'] == 'SSB' else 'dead_catch_tonnes'
        actual = float(projections[key][field]);expected = float(r['published_bootstrap_median_tonnes'])
        comparisons.append(dict(scenario=key[0], year=key[1], quantity=r['quantity'], quadra_deterministic_tonnes=actual, published_bootstrap_median_tonnes=expected, difference_tonnes=actual-expected, relative_difference=(actual-expected)/expected))
    with (OUT/'published_projection_comparison.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=comparisons[0]);writer.writeheader();writer.writerows(comparisons)
    sensitivity_comparisons = []
    for r in read('projection_sensitivities.csv'):
        base = projections[int(r['baseline_scenario']), int(r['year'])]
        sensitivity_comparisons.append(dict(sensitivity=r['sensitivity'], baseline_scenario=int(r['baseline_scenario']), year=int(r['year']), spawning_biomass_difference_tonnes=float(r['spawning_biomass_tonnes'])-float(base['spawning_biomass_tonnes']), dead_catch_difference_tonnes=float(r['dead_catch_tonnes'])-float(base['dead_catch_tonnes'])))
    if len(sensitivity_comparisons) != 30:
        raise ValueError('Expected three ten-year source interpretation sensitivities')
    with (OUT/'projection_sensitivity_comparison.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=sensitivity_comparisons[0]);writer.writeheader();writer.writerows(sensitivity_comparisons)
    (OUT/'projection_comparison_summary.json').write_text(json.dumps(dict(status='deterministic reconstruction; published bootstrap medians are not validation targets', scenarios=5, years=[2022,2031], quantities_compared=100, maximum_absolute_relative_difference=max(abs(r['relative_difference']) for r in comparisons), unresolved=['Original scenario forecast files and 100 bootstrap fitted models were not supplied.','Methods specify 2008-2010 reference SPR; table labels 1998-2000.','Status quo uses matching equilibrium annual F with fixed apical rates; realized pooled annual F varies with age composition.','Scenario 1 uses the literal dynamic baseline; equilibrium 20% B0 is a different target.']),indent=2)+'\n')
    if not dynamic['passed'] or not spr_result['passed']:
        raise SystemExit('Dynamic unfished comparison failed')
    print(f"All {len(checks)} dynamic unfished comparisons pass. All 47 historical SPR checks pass. Five projection scenarios written; published bootstrap differences retained.")
if __name__ == '__main__':
    main()
