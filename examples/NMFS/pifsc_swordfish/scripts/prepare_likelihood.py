#!/usr/bin/env python3
"""Prepare likelihood comparison targets, retaining original observations as model inputs."""
import argparse
from pathlib import Path
from prepare_model import ROOT,save_if_changed,report_table
from import_ss3 import read_lines,extract
from swordfish_inputs import load_inputs


def active(path):
    """Return nonempty source lines after removing comments to check supported starter and forecast settings.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.

    Returns
    -------
    Nonempty source lines with comments removed.

    Notes
    -----
    The target is the joint penalized likelihood with configured phase lambdas, profiled Q and supplied composition weighting. Forecast settings are checked for the supported disabled-projection input, not inferred from narrative scenario descriptions.
    """
    return [s.split('#',1)[0].split() for s in read_lines(path) if s.split('#',1)[0].strip()]


def prepare(raw,output):
    """Validate the supported starter/forecast likelihood settings and generate named likelihood/index/length comparison targets.

    Parameters
    ----------
    raw
        Directory containing the eight original assessment inputs and reference report.
    output
        Destination for generated artifacts; the function summary specifies a file or directory.

    Notes
    -----
    The target is the joint penalized likelihood with configured phase lambdas, profiled Q and supplied composition weighting. Forecast settings are checked for the supported disabled-projection input, not inferred from narrative scenario descriptions.
    """
    model=load_inputs(raw)
    # The implemented objective is scoped to these supplied starter/forecast options.
    starter=active(raw/'starter.ss');forecast=active(raw/'forecast.ss')
    if starter[9]!=['1']:raise ValueError('Only supplied enabled selectivity soft bounds supported')
    if forecast[6]!=['0'] or forecast[7]!=['0']:raise ValueError('Likelihood comparison excludes active projections')
    if model['control']['ballpark_year']>=0:raise ValueError('Active F-ballpark penalty unsupported')
    ls=read_lines(raw/'Report.sso')
    _,rows=extract(ls,'LIKELIHOOD','Component logL*Lambda',r'^\S+\s+[-+\d]')
    components=[(r[0],r[1]) for r in rows]
    indices=report_table(ls,'INDEX_2','Fleet Fleet_name')
    lengths=report_table(ls,'FIT_LEN_COMPS','Fleet Fleet_Name')
    lines=['SWORDFISH_LIKELIHOOD_TARGETS_V1',str(len(components))]
    lines += [f'{name} {v}' for name,v in components]
    lines += [str(len(indices))]
    for r in indices:
        lines.append(' '.join(r[k] for k in ['Yr','Seas','Fleet','Vuln_bio','Exp','Calc_Q','SE','Dev','Like+log(s)']))
    lines += [str(len(lengths))]
    for r in lengths:
        lines.append(' '.join(r[k] for k in ['Yr','Seas','Fleet','Nsamp_adj','Like','All_obs_mean','All_exp_mean']))
    save_if_changed(output/'likelihood_targets.txt','\n'.join(lines)+'\n')
    print(f'Prepared {len(components)} objective components, {len(indices)} index fits and {len(lengths)} length fits as targets only.')


def main():
    """Execute the prepare likelihood command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    The target is the joint penalized likelihood with configured phase lambdas, profiled Q and supplied composition weighting. Forecast settings are checked for the supported disabled-projection input, not inferred from narrative scenario descriptions.
    """
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,default=ROOT/'data/raw');p.add_argument('--output',type=Path,default=ROOT/'build')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True);prepare(args.raw,args.output)

if __name__=='__main__':main()
