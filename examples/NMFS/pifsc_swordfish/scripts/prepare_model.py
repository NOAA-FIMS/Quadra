#!/usr/bin/env python3
"""Export original-input schema and generate C++ data + conditional SS3 replay fixtures."""
import argparse
import json
from pathlib import Path
from swordfish_inputs import load_inputs
from import_ss3 import read_lines, section, extract

ROOT = Path(__file__).resolve().parents[1]


def save_if_changed(path, text):
    # Avoid recompiling the generated C++ factory when inputs are unchanged.
    """Write a generated artifact only when content differs, preserving header timestamps to avoid unnecessary recompilation.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.
    text
        Source/document text; output functions preserve intended line breaks.

    Notes
    -----
    Generated C++ factories contain original observations and controls. Report trajectories are comparison fixtures, not drivers of the independent population dynamics. Do not hand-edit generated build headers; regenerate from the validated sources.
    """
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


def cpp(x):
    """Serialize supported scalar, string and nested-list input values as C++ initializer expressions.

    Parameters
    ----------
    x
        Scalar/string/list value recursively converted to a C++ initializer.

    Returns
    -------
    C++ initializer text, recursively serialized without executing source comments.

    Notes
    -----
    Generated C++ factories contain original observations and controls. Report trajectories are comparison fixtures, not drivers of the independent population dynamics. Do not hand-edit generated build headers; regenerate from the validated sources.
    """
    if isinstance(x, str): return json.dumps(x)
    if isinstance(x, bool): return 'true' if x else 'false'
    if isinstance(x, (int,float)): return repr(x)
    return '{' + ','.join(cpp(v) for v in x) + '}'


def generate_header(model, output, include_observations=True, observation_kinds=None, factory_name=None):
    """Generate a C++ InputData factory from parsed original observations and controls, with optional observation subsets for audits.

    Parameters
    ----------
    model
        Validated Python dictionary containing original data, controls and provenance.
    output
        Destination for generated artifacts; the function summary specifies a file or directory.
    include_observations
        Whether to emit observation arrays rather than only biological controls.
    observation_kinds
        Optional observation-family filter used for audit-specific factories.
    factory_name
        Optional generated C++ factory name; omitted uses the supported default.

    Notes
    -----
    Generated C++ factories contain original observations and controls. Report trajectories are comparison fixtures, not drivers of the independent population dynamics. Do not hand-edit generated build headers; regenerate from the validated sources.
    """
    data, ctl = model['data'], model['control']
    d = data['dimensions']
    p = {r['name']:r['columns'][2] for r in ctl['parameters']}
    lines = ['// Generated from original .dat/.ctl by prepare_model.py. Do not edit.',
             '#pragma once', '#include "swordfish/data.hpp"',
             'inline swordfish::InputData imported_swordfish_data() {',
             'swordfish::InputData d;']
    def assign(name,value):
        """Append one named field assignment with a typed serialized initializer to the generated C++ factory.

        Parameters
        ----------
        name
            Exact section, output, field or stage identifier; matching is case sensitive.
        value
            Value to serialize or include in the generated fixture record.
        """
        lines.append(f'd.{name} = {cpp(value)};')
    assign('dimensions',[d['start_year'],d['end_year'],d['seasons'],d['subseasons'],
                         d['plus_group_age'],d['areas'],d['spawn_month'],d['season_months']])
    assign('biology.natural_mortality',[v for sex in ctl['natural_mortality'] for v in sex])
    for sex,suffix in enumerate(['Fem','Mal']):
        names = ['L_at_Amin','L_at_Amax','VonBert_K','CV_young','CV_old','Wtlen_1','Wtlen_2']
        vals = [p[f'{name}_{suffix}_GP_1'] for name in names]
        assign(f'biology.growth[{sex}]',ctl['growth_reference_ages'] + vals)
    for key,name in [('maturity_l50','Mat50%_Fem_GP_1'),('maturity_slope','Mat_slope_Fem_GP_1'),
                     ('female_fraction','FracFemale_GP_1')]:
        assign('biology.'+key,p[name])
    assign('biology.first_mature_age',ctl['first_mature_age'])
    assign('biology.settlement_month',ctl['settlement'][1])
    assign('biology.settlement_age',ctl['settlement'][3])
    for f in data['fleets']:
        lines.append('d.fleets.push_back('+cpp([f[k] for k in
                     ['id','type','timing','area','catch_units','catch_multiplier','name']])+');')
    for name,source,keys in [
        ('catches','catches',['year','season','fleet','value','log_se']),
        ('indices','indices',['year','month','fleet','value','log_se','use']),
        ('compositions','length_compositions',['year','month','fleet','sex','partition','sample_size','use','values'])] if include_observations else []:
        if observation_kinds is None or name in observation_kinds:
            assign(name,[[row[k] for k in keys] for row in data[source]])
    for name in ['population_length_bins','length_bins','index_specs','length_controls']:
        assign(name,data[name])
    for name in ['blocks','q_setup','size_patterns','age_patterns','variance_adjustments','lambda_changes']:
        assign(name,ctl[name])
    for row in ctl['parameters']:
        lines.append('d.parameters.push_back('+cpp([row['name'],row['columns'],row['source_line']])+');')
    for k in ['start','end','phase']: assign('recdev_'+k,ctl['recruitment'][k])
    assign('recruitment_advanced',ctl['recruitment']['advanced'])
    for name in ['f_method','hybrid_iterations','max_lambda_phase','sd_offset','max_f','plus_group_growth']:
        assign(name,ctl[name])
    lines += ['d.validate(); return d;', '}']
    text='\n'.join(lines)+'\n'
    if not include_observations:
        text=text.replace('imported_swordfish_data()', 'imported_swordfish_biology_data()')
    if factory_name is not None:
        text=text.replace('imported_swordfish_data()',factory_name+'()')
    save_if_changed(output,text)


def report_table(lines,name,prefix):
    """Locate a version-specific rectangular Report.sso section and return its header and token rows for comparison fixtures.

    Parameters
    ----------
    lines
        Source report lines after the documented encoding/newline normalization.
    name
        Exact section, output, field or stage identifier; matching is case sensitive.
    prefix
        Exact expected report header prefix.

    Returns
    -------
    Rectangular report rows keyed by the checked header.

    Notes
    -----
    Generated C++ factories contain original observations and controls. Report trajectories are comparison fixtures, not drivers of the independent population dynamics. Do not hand-edit generated build headers; regenerate from the validated sources.
    """
    block = section(lines,name)
    h = next(i for i,s in enumerate(block) if s.lstrip().startswith(prefix))
    header = block[h].split()
    rows=[]
    for s in block[h+1:]:
        t=s.split()
        if not t or t[0].startswith('#'): break
        if len(t) != len(header):
            raise ValueError(f'{name}: unsupported row width {len(t)} vs {len(header)}')
        rows.append(dict(zip(header,t)))
    if not rows: raise ValueError(f'{name}: empty report table')
    return rows


def generate_replay(model, raw, output):
    """Report states/rates are fixtures only, never model observation inputs.

    Parameters
    ----------
    model
        Validated Python dictionary containing original data, controls and provenance.
    raw
        Directory containing the eight original assessment inputs and reference report.
    output
        Destination for generated artifacts; the function summary specifies a file or directory.

    Notes
    -----
    Generated C++ factories contain original observations and controls. Report trajectories are comparison fixtures, not drivers of the independent population dynamics. Do not hand-edit generated build headers; regenerate from the validated sources.
    """
    ls=read_lines(raw/'Report.sso')
    naa=report_table(ls,'NUMBERS_AT_AGE','Area Bio_Pattern')
    faa=report_table(ls,'F_AT_AGE','Area Fleet Sex')
    # CATCH_AT_AGE repeats the unused XX header: index numeric columns explicitly.
    block=section(ls,'CATCH_AT_AGE')
    catch={}
    for line in block[1:]:
        t=line.split()
        if not t: break
        if len(t)!=27: raise ValueError('Unexpected catch-at-age width')
        if t[10]=='TIME' and t[5]=='dead':
            key=(int(t[7]),int(t[8]),int(t[2]),int(t[1]))
            if key in catch: raise ValueError('Duplicate catch-at-age fixture')
            catch[key]=list(map(float,t[11:]))
    rates={}
    for row in faa:
        if row['Era']!='TIME': continue
        key=tuple(int(row[k]) for k in ['Yr','Seas','Sex','Fleet'])
        if key in rates: raise ValueError('Duplicate F fixture')
        rates[key]=[float(row[str(a)]) for a in range(16)]
    states={}
    for row in naa:
        if row['Era']!='TIME': continue
        key=(int(row['Yr']),int(row['Seas']),int(row['Sex']),row['Beg/Mid'])
        if key in states: raise ValueError('Duplicate population fixture')
        states[key]=[float(row[str(a)]) for a in range(16)]
    wt=report_table(ls,'MEAN_BODY_WT(Begin)','Morph Yr Seas')
    weights={(int(r['Yr']),int(r['Seas']),int(r['Morph'])):[float(r[str(a)]) for a in range(16)] for r in wt}
    bio=report_table(ls,'Biology_at_age_in_endyr','Seas Morph')
    mature_weights={int(r['int_Age']):float(r['Mat*Fecund']) for r in bio if r['Seas']=='2' and r['Sex']=='1'}
    # Preserve runtime biology separately from original controls.
    growth=report_table(ls,'Growth_Parameters','Count Yr Sex')
    header, rows = extract(ls,'DERIVED_QUANTITIES','Label Value',r'^\S+\s+[-+\d]')
    derived = {r[0]:float(r[1]) for r in rows}
    time={}
    for r in report_table(ls,'TIME_SERIES','Area Yr Era Seas'):
        if r['Era']=='TIME': time[(int(r['Yr']),int(r['Seas']))]=r
    d=model['data']['dimensions']
    values=['SWORDFISH_CONDITIONAL_REPLAY_V2',str(d['start_year']),str(d['end_year'])]
    values += [str(mature_weights[a]) for a in range(16)]
    # Year/season order followed by state, midstate, body weights, fleet rates,
    # catches, then recruitment, biomass and SSB targets. All rates are from SS3.
    for y in range(d['start_year'],d['end_year']+1):
        for s in range(1,5):
            values += [str(y),str(s)]
            for sex in (1,2):
                values += list(map(str,states[y,s,sex,'B']))
                values += list(map(str,states[y,s,sex,'M']))
                values += list(map(str,weights[y,s,sex]))
                for f in range(1,20):
                    values += list(map(str,rates[y,s,sex,f]))
                    values += list(map(str,catch[y,s,sex,f]))
            row=time[y,s]
            values += [str(derived[f'Recr_{y}']),row['Bio_all'],row['Bio_smry'],
                       row['SpawnBio'] if row['SpawnBio']!='_' else '-1']
            values += [row[f'F:_{f}'] for f in range(1,20)]
    values += [str(len(bio))]
    for r in bio:
        values += [r[k] for k in ['Sex','Seas','int_Age','Len_Beg','Len_Mid']]
    values += [str(len(growth))]
    for r in growth:
        values += [r[k] for k in ['Sex','A1','A2','L_a_A1','L_a_A2','K','CVmin','CVmax','WtLen1','WtLen2']]
    save_if_changed(output,'\n'.join(values)+'\n')
    header, rows = extract(ls,'PARAMETERS','Num Label Value',r'^\d+\s')
    report_fixed = {r[1]:float(r[2]) for r in rows if float(r[4])<0}
    differences=[]
    for p in model['control']['parameters']:
        if p['name'] in report_fixed and p['columns'][6]<0:
            v=report_fixed[p['name']]
            if v!=p['columns'][2]:
                differences.append(dict(parameter=p['name'],control_initial=p['columns'][2],report_value=v))
    return differences


def main():
    """Execute the prepare model command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    Generated C++ factories contain original observations and controls. Report trajectories are comparison fixtures, not drivers of the independent population dynamics. Do not hand-edit generated build headers; regenerate from the validated sources.
    """
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--raw',type=Path,default=ROOT/'data/raw')
    ap.add_argument('--output',type=Path,default=ROOT/'build')
    args=ap.parse_args()
    model=load_inputs(args.raw)
    args.output.mkdir(parents=True,exist_ok=True)
    save_if_changed(args.output/'swordfish_inputs.json',json.dumps(model,indent=2)+'\n')
    generate_header(model,args.output/'swordfish_input.hpp')
    differences=generate_replay(model,args.raw,args.output/'ss3_replay.txt')
    summary=dict(raw_catch_rows=len(model['data']['catches']),raw_index_rows=len(model['data']['indices']),
                 raw_length_samples=len(model['data']['length_compositions']),
                 included_length_samples=sum(x['use'] for x in model['data']['length_compositions']),
                 control_parameters=len(model['control']['parameters']),
                 fixed_parameter_differences=differences,
                 status='Original inputs parsed; report fixtures enable conditional dynamics validation only')
    save_if_changed(args.output/'preparation_summary.json',json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
