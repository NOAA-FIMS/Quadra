#!/usr/bin/env python3
"""Prepare independent-population inputs and SS3 validation targets."""
import argparse
from pathlib import Path
from prepare_model import ROOT,generate_header,save_if_changed,report_table
from swordfish_inputs import load_inputs
from import_ss3 import read_lines,extract


def prepare(raw,output):
    """Generate original-catch population inputs, early/main fitted deviation mappings and initialization targets from the supplied report.

    Parameters
    ----------
    raw
        Directory containing the eight original assessment inputs and reference report.
    output
        Destination for generated artifacts; the function summary specifies a file or directory.

    Notes
    -----
    The early age deviations are mapped to their recruitment years. The final main deviation is dependent under the zero-sum constraint in fitting; report values are used only to construct validation fixtures and starting values.
    """
    model=load_inputs(raw)
    generate_header(model,output/'swordfish_population_input.hpp',observation_kinds={'catches'},
                    factory_name='imported_swordfish_population_data')
    # This factory contains original catch observations, no report-derived states.
    ls=read_lines(raw/'Report.sso')
    _,rows=extract(ls,'PARAMETERS','Num Label Value',r'^\d+\s')
    first=model['data']['dimensions']['start_year']
    last=model['data']['dimensions']['end_year']
    deviations={}
    for row in rows:
        label=row[1]
        if label.startswith('Early_InitAge_'):year=first-int(label.split('_')[-1])
        elif label.startswith(('Early_RecrDev_','Main_RecrDev_')):year=int(label.split('_')[-1])
        else:continue
        if year in deviations:raise ValueError('Duplicate recruitment deviation')
        deviations[year]=row[2]
    expected=set(range(int(model['control']['recruitment']['advanced'][0]),last+1))
    if set(deviations)!=expected:raise ValueError('Incomplete fitted recruitment deviation set')
    save_if_changed(output/'population_fitted_deviations.txt',
                    '\n'.join(['SWORDFISH_DEVIATIONS_V1',str(len(deviations))]+[f'{y} {deviations[y]}' for y in sorted(deviations)])+'\n')
    _,rows=extract(ls,'DERIVED_QUANTITIES','Label Value',r'^\S+\s+[-+\d]')
    derived={r[0]:r[1] for r in rows}
    values=['SWORDFISH_INITIAL_TARGETS_V1']
    values += [derived[k] for k in ['Recr_Virgin','SSB_Virgin','Recr_Initial','SSB_Initial']]
    naa=report_table(ls,'NUMBERS_AT_AGE','Area Bio_Pattern')
    for era,year in [('VIRG',first-2),('INIT',first-1),('TIME',first)]:
        for sex in (1,2):
            r=next(r for r in naa if r['Era']==era and r['Yr']==str(year) and r['Sex']==str(sex) and r['Seas']=='1' and r['Beg/Mid']=='B')
            values += [r[str(a)] for a in range(16)]
    save_if_changed(output/'initialization_targets.txt','\n'.join(values)+'\n')
    print(f'Prepared original catch inputs and {len(deviations)} fitted recruitment deviations; no population states used as inputs.')


def main():
    """Execute the prepare population command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    The early age deviations are mapped to their recruitment years. The final main deviation is dependent under the zero-sum constraint in fitting; report values are used only to construct validation fixtures and starting values.
    """
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,default=ROOT/'data/raw');p.add_argument('--output',type=Path,default=ROOT/'build')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True);prepare(args.raw,args.output)

if __name__=='__main__':main()
