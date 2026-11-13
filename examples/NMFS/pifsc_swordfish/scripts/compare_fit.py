#!/usr/bin/env python3
"""Compare saved Quadra fitted annual outputs to the supplied SS3 report."""
import csv
import json
from pathlib import Path
from prepare_model import ROOT,report_table
from import_ss3 import read_lines,extract


def main():
    """Execute the compare fit command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    Compare independent fitted outputs against rounded SS3 report quantities using quantity-specific absolute and relative tolerances. A nonzero exit means at least one required comparison failed; passing does not imply all SS3 configurations are supported.
    """
    build=ROOT/'build';lines=read_lines(ROOT/'data/raw/Report.sso')
    _,rows=extract(lines,'DERIVED_QUANTITIES','Label Value',r'^\S+\s+[-+\d]')
    derived={r[0]:float(r[1]) for r in rows}
    time={int(r['Yr']):r for r in report_table(lines,'TIME_SERIES','Area Yr Era Seas') if r['Era']=='TIME' and r['Seas']=='1'}
    spr={int(r['Yr']):r for r in report_table(lines,'SPR_SERIES','Yr Era Bio_all') if r['Era']=='TIME'}
    annual=list(csv.DictReader((build/'quadra_fitted_annual.csv').open()))
    if {int(r['year']) for r in annual}!=set(time) or len(annual)!=len(time):raise ValueError('Incomplete fitted annual series')
    summary={};output=[]
    for row in annual:
        y=int(row['year'])
        targets={'spawning_biomass':derived[f'SSB_{y}'],'recruitment_thousands':derived[f'Recr_{y}'],
                 'age1plus_biomass_january':float(time[y]['Bio_smry']),'annual_F_SS3_units4_ages1_10':float(spr[y]['F=Z-M'])}
        result={'year':y}
        for quantity,target in targets.items():
            fitted=float(row[quantity]);error=abs(fitted-target)
            atol=1e-6 if quantity.startswith('annual_F') else (0.002 if quantity.startswith('recruitment') else 0.7)
            passed=error<=atol+2e-5*abs(target)
            result.update({quantity+'_quadra':fitted,quantity+'_ss3':target,quantity+'_difference':fitted-target})
            metric=summary.setdefault(quantity,{'comparisons':0,'max_absolute_error':0,'passed':True})
            metric['comparisons']+=1;metric['max_absolute_error']=max(metric['max_absolute_error'],error);metric['passed'] &= passed
        output.append(result)
    with (build/'quadra_fit_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    (build/'quadra_fit_comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not all(r['passed'] for r in summary.values()):raise SystemExit(1)

if __name__=='__main__':main()
