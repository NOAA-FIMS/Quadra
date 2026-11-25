"""Independent native-schema mapping for conversion parity tests."""
from native_input import DIMENSIONS, GROWTH, RECRUITMENT, UNITS, encode


def convert(source, values):
    """Map parsed observations/settings and explicit native starts into the named schema.

    Existing numeric settings are preserved, including selectivity block
    functions, bias ramp, observation inclusion and likelihood weights. Native
    starts may come from a complete Quadra parameter CSV. No fitted trajectory
    or SS3 report quantities become observations or population drivers.
    """
    d,c=source['data'],source['control'];p={r['name']:r['columns'][2] for r in c['parameters']}
    names=set(p)|{f'RecruitmentDev_{y}' for y in range(int(c['recruitment']['advanced'][0]),c['recruitment']['end']+1)}
    if values is not None and set(values)!=names:raise ValueError('Initial parameter CSV must be complete and contain no unknown names')
    if values is not None:p.update({k:values[k] for k in p})
    dims=d['dimensions']
    dim_values=[dims[k] for k in ['start_year','end_year','seasons','subseasons','plus_group_age','areas','spawn_month','season_months']]
    biology=dict(natural_mortality=[v for sex in c['natural_mortality'] for v in sex],growth=[],
                 maturity_l50=p['Mat50%_Fem_GP_1'],maturity_slope=p['Mat_slope_Fem_GP_1'],female_fraction=p['FracFemale_GP_1'],
                 first_mature_age=c['first_mature_age'],settlement_month=c['settlement'][1],settlement_age=c['settlement'][3])
    for suffix in ['Fem','Mal']:
        vals=c['growth_reference_ages']+[p[f'{name}_{suffix}_GP_1'] for name in ['L_at_Amin','L_at_Amax','VonBert_K','CV_young','CV_old','Wtlen_1','Wtlen_2']]
        biology['growth'].append(dict(zip([k for k,_ in GROWTH],vals)))
    model=dict(dimensions=dict(zip([k for k,_ in DIMENSIONS],dim_values)),biology=biology)
    for key in ['fleets','catches','indices']:
        keys={'fleets':['id','type','timing','area','catch_units','catch_multiplier','name'],
              'catches':['year','season','fleet','value','log_se'],
              'indices':['year','month','fleet','value','log_se','use']}[key]
        model[key]=[{k:r[k] for k in keys} for r in d[key]]
    keys=['year','month','fleet','sex','partition','sample_size','use','values']
    model['compositions']=[{k:r[k] for k in keys} for r in d['length_compositions']]
    for key in ['population_length_bins','length_bins','index_specs','length_controls']:model[key]=d[key]
    for key in ['blocks','q_setup','size_patterns','age_patterns','variance_adjustments','lambda_changes','f_method','hybrid_iterations','max_lambda_phase','sd_offset','max_f','plus_group_growth']:model[key]=c[key]
    rec=c['recruitment'];rec_values=[rec[k] for k in ['start','end','phase']]+rec['advanced']
    model['recruitment']={k:int(v) if kind=='i' else v for (k,kind),v in zip(RECRUITMENT,rec_values)}
    parameters=[]
    for spec in c['parameters']:
        col=spec['columns'];parameters.append(dict(name=spec['name'],lower=col[0],upper=col[1],value=p[spec['name']],phase=int(col[6]),prior_type=int(col[5]),block=int(col[12]) if len(col)>12 else 0,block_function=int(col[13]) if len(col)>13 else 0))
    dev=[dict(year=y,value=values[f'RecruitmentDev_{y}'] if values is not None else 0.) for y in range(int(rec['advanced'][0]),rec['end']+1)]
    result=dict(schema='quadra.swordfish.v1',units=UNITS,model=model,parameters=parameters,recruitment_deviations=dev)
    encode(result)
    return result
