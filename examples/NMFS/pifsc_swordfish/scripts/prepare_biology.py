#!/usr/bin/env python3
"""Prepare report-only parameter values and biological prediction fixtures."""
import argparse
import csv
from pathlib import Path
import re
from import_ss3 import read_lines, section, extract
from swordfish_inputs import load_inputs
from prepare_model import ROOT, generate_header, report_table, save_if_changed


def prepare(raw, output):
    """Generate biological input/parameter mappings and report fixtures, validating canonical parameter order and reporting the fixed male-growth discrepancy.

    Parameters
    ----------
    raw
        Directory containing the eight original assessment inputs and reference report.
    output
        Destination for generated artifacts; the function summary specifies a file or directory.

    Notes
    -----
    Parameter mapping is checked by canonical name/order and bounds. The fixed male growth length differs between control initialization (83.2) and fitted report (83.4); the fitted-value fixture preserves the report value explicitly.
    """
    model=load_inputs(raw)
    generate_header(model,output/'swordfish_biology_input.hpp',include_observations=False)
    ls=read_lines(raw/'Report.sso')
    _, rows=extract(ls,'PARAMETERS','Num Label Value',r'^\d+\s')
    fitted=[r for r in rows if not r[1].startswith(('Early_','Main_','ForeRecr_'))]
    specs=model['control']['parameters']
    if len(fitted)!=len(specs): raise ValueError('Report/input parameter count mismatch')
    mapping=[]
    for spec,row in zip(specs,fitted):
        name=spec['name']; label=row[1]
        if name.startswith('SizeSel_fleet_'):
            f=int(name.split('_')[2])
            if f'({f})' not in label: raise ValueError(f'Selectivity fleet order mismatch: {name} / {label}')
        elif name.startswith('LnQ_fleet_'):
            if f"({name.split('_')[-1]})" not in label: raise ValueError('Q fleet order mismatch')
        elif not label.startswith(name): raise ValueError(f'Parameter order mismatch: {name} / {label}')
        if float(row[5])!=spec['columns'][0] or float(row[6])!=spec['columns'][1]:
            raise ValueError(f'Parameter bounds mismatch: {name}')
        mapping.append((name,label,row[2]))
    parameter_lines=['SWORDFISH_FITTED_PARAMETERS_V1',str(len(mapping))]
    parameter_lines += [f'{name} {v}' for name,_,v in mapping]
    save_if_changed(output/'biology_fitted_parameters.txt','\n'.join(parameter_lines)+'\n')
    values=['SWORDFISH_BIOLOGY_FIXTURE_V1',str(len(mapping))]
    for name,_,v in mapping: values += [name,v]
    cases=[]
    def case(kind,year,season,sub,sex,fleet,age,bin,value):
        """Append one labeled biological fixture row with year, timing, sex, fleet, age, bin and reference value.

        Parameters
        ----------
        kind
            Comparison-fixture quantity label.
        year
            Calendar year of the fixture.
        season
            One-based model quarter.
        sub
            Source subseason code.
        sex
            Source fixture sex code (1 female, 2 male), not C++ zero-based sex.
        fleet
            One-based source fleet identifier, or documented fixture sentinel.
        age
            Model age, with age 15 the plus group.
        bin
            Source length-bin index or a documented non-bin sentinel.
        value
            Value to serialize or include in the generated fixture record.
        """
        cases.append(' '.join(map(str,[kind,year,season,sub,sex,fleet,age,bin,value])))
    # Report biological ages and derived weights, including the plus group.
    for r in report_table(ls,'Biology_at_age_in_endyr','Seas Morph'):
        for kind,col in [('length_begin','Len_Beg'),('length_mid','Len_Mid'),
                         ('sd_begin','SD_Beg'),('sd_mid','SD_Mid'),
                         ('weight_begin','Wt_Beg'),('weight_mid','Wt_Mid')]:
            case(kind,2021,r['Seas'],0,r['Sex'],0,r['int_Age'],0,r[col])
        if r['Sex']=='1':case('mature_weight',2021,r['Seas'],0,1,0,r['int_Age'],0,r['Mat*Fecund'])
        for f in range(1,28):
            case('selected_weight',2021,r['Seas'],0,r['Sex'],f,r['int_Age'],0,r[f'SelWt:_{f}'])
            # Pattern 0 reports the continuous mean rather than bin midpoints.
            if model['control']['size_patterns'][f-1][0]!=0:
                case('selected_length',2021,r['Seas'],0,r['Sex'],f,r['int_Age'],0,r[f'Len:_{f}'])
    # ALK tables are transposed: descending length edges, 16 age columns.
    block=section(ls,'AGE_LENGTH_KEY')
    current=None; count=0
    for line in block:
        match=re.match(r'\s*Seas:\s*(\d+) Sub_Seas:\s*(\d+)\s+Morph:\s*(\d+)',line)
        if match:
            if current is not None and count!=51: raise ValueError('Incomplete ALK table')
            current=tuple(map(int,match.groups()));count=0;continue
        t=line.split()
        if current and len(t)==17 and re.fullmatch(r'\d+(?:\.\d+)?',t[0]):
            season,sub,sex=current; b=int((float(t[0])-10)/5)
            for a,v in enumerate(t[1:]):case('alk',2021,season,sub,sex,0,a,b,v)
            count+=1
    if count!=51:raise ValueError('Incomplete final ALK table')
    # Explicit year-level curves check time blocks and mirror propagation.
    for r in report_table(ls,'LEN_SELEX','Factor Fleet Yr Sex'):
        if r['Factor']=='Lsel' and 1975<=int(r['Yr'])<=2021:
            for b in range(51):
                case('length_selectivity',r['Yr'],0,0,r['Sex'],r['Fleet'],0,b,r[str(12.5+5*b)])
    for r in report_table(ls,'AGE_SELEX','Factor Fleet Yr Seas'):
        if r['Factor']=='Asel2' and 1975<=int(r['Yr'])<=2021:
            for a in range(16):case('age_selectivity',r['Yr'],r['Seas'],0,r['Sex'],r['Fleet'],a,0,r[str(a)])
    # Compare effective parameters for every historical year, not just plotted curves.
    for line in section(ls,'selparm(Size)_By_Year_after_adjustments')[1:]:
        t=line.split()
        if not t:break
        f,y=int(t[0]),int(t[1])
        if 1975<=y<=2021:
            for i,v in enumerate(t[3:],1):case('effective_parameter',y,0,0,1,f,0,i,v)
    values += [str(len(cases)),*cases]
    save_if_changed(output/'ss3_biology_fixtures.txt','\n'.join(values)+'\n')
    with (output/'biology_parameter_mapping.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['canonical_name','report_label','report_value']);w.writerows(mapping)
    print(f'Prepared {len(cases)} biological comparisons using {len(mapping)} explicitly mapped fitted parameters.')


def main():
    """Execute the prepare biology command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    Parameter mapping is checked by canonical name/order and bounds. The fixed male growth length differs between control initialization (83.2) and fitted report (83.4); the fitted-value fixture preserves the report value explicitly.
    """
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,default=ROOT/'data/raw')
    p.add_argument('--output',type=Path,default=ROOT/'build')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    prepare(args.raw,args.output)

if __name__=='__main__':main()
