"""Shared schema checks for research/release tools and independent tests.

Production JSON intake and validation live in C++; this module is not a runner.
"""
import json
import math

DIMENSIONS = [('first_year','i'),('last_year','i'),('seasons','i'),('subseasons','i'),('plus_age','i'),('areas','i'),('spawn_month','d'),('season_months',('a',4,'d'))]
GROWTH = [(k,'d') for k in ['age1','age2','length1','length2','k','cv_young','cv_old','weight_coefficient','weight_power']]
BIOLOGY = [('natural_mortality',('a',32,'d')),('maturity_l50','d'),('maturity_slope','d'),('female_fraction','d'),('first_mature_age','i'),('settlement_month','d'),('settlement_age','d')]
RECORDS = {
    'fleets': [(k,'i') for k in ['id','type','timing','area','catch_units','catch_multiplier']]+[('name','s')],
    'catches': [(k,'i') for k in ['year','season','fleet']]+[('value','d'),('log_se','d')],
    'indices': [('year','i'),('month','d'),('fleet','i'),('value','d'),('log_se','d'),('use','b')],
    'compositions': [('year','i'),('month','d'),('fleet','i'),('sex','i'),('partition','i'),('sample_size','d'),('use','b'),('values',('a',102,'d'))],
}
TABLES = {'index_specs':(4,'i'),'size_patterns':(4,'i'),'age_patterns':(4,'i'),
          'length_controls':(7,'d'),'q_setup':(6,'d'),'variance_adjustments':(3,'d'),'lambda_changes':(5,'d')}
SCALARS = [(k,'i') for k in ['f_method','hybrid_iterations','max_lambda_phase','sd_offset']]+[('max_f','d'),('plus_group_growth','d')]
RECRUITMENT = [(k,'i') for k in ['main_start','main_end','phase','early_start','early_phase','forecast_phase']]+[(k,'d') for k in ['forecast_lambda','last_early_no_bias','first_full_bias','last_full_bias','first_recent_no_bias','maximum_bias']]+[('cycles','i'),('minimum_deviation','d'),('maximum_deviation','d'),('read_deviations','i')]
PARAMETER = [('name','s'),('lower','d'),('upper','d'),('value','d'),('phase','i'),('prior_type','i'),('block','i'),('block_function','i')]
UNITS = dict(numbers='thousands of fish', weight='kg', biomass='metric tonnes', length='cm EFL', mortality='annual instantaneous rate')
BIO_NAMES = ['L_at_Amin_Fem_GP_1','L_at_Amax_Fem_GP_1','VonBert_K_Fem_GP_1','CV_young_Fem_GP_1','CV_old_Fem_GP_1','Wtlen_1_Fem_GP_1','Wtlen_2_Fem_GP_1','Mat50%_Fem_GP_1','Mat_slope_Fem_GP_1','Eggs/kg_inter_Fem_GP_1','Eggs/kg_slope_wt_Fem_GP_1','L_at_Amin_Mal_GP_1','L_at_Amax_Mal_GP_1','VonBert_K_Mal_GP_1','CV_young_Mal_GP_1','CV_old_Mal_GP_1','Wtlen_1_Mal_GP_1','Wtlen_2_Mal_GP_1','CohortGrowDev','FracFemale_GP_1']
PARAMETER_NAMES = set(BIO_NAMES+['SR_LN(R0)','SR_BH_steep','SR_sigmaR','SR_regime','SR_autocorr','InitF_seas_1_flt_6']+[f'LnQ_fleet_{f}' for f in range(20,28)]+[f'SizeSel_fleet_{f}_p{i}' for f,n in [(1,6),(2,2),(3,6),(4,2),(6,6),(7,6),(8,6)] for i in range(1,n+1)]+[f'SizeSel_fleet_1_p{i}_block_2_1' for i in range(1,5)]+[f'SizeSel_fleet_2_p{i}_block_1_1' for i in range(1,3)])


def exact_keys(value, keys, label):
    """Reject missing/unknown fields instead of silently accepting misspelled settings."""
    if not isinstance(value, dict) or set(value) != set(keys):
        raise ValueError(f'{label}: expected exactly {sorted(keys)}')


def atom(value, kind):
    """Validate a scalar/fixed array and encode finite typed values for the C++ reader."""
    if isinstance(kind, tuple):
        if not isinstance(value, list) or len(value) != kind[1]:
            raise ValueError(f'Expected array of length {kind[1]}')
        return ' '.join(atom(v,kind[2]) for v in value)
    if kind == 's':
        if not isinstance(value,str) or any(ord(c)<32 for c in value):
            raise ValueError('Expected string without control characters')
        return json.dumps(value, ensure_ascii=False)
    if kind == 'b':
        if type(value) is not bool: raise ValueError('Expected boolean')
        return str(int(value))
    if type(value) not in (int,float) or not math.isfinite(value):
        raise ValueError('Expected finite number')
    if kind == 'i':
        if type(value) is not int: raise ValueError('Expected integer')
        if not -(2**31) <= value < 2**31: raise ValueError('Integer out of range')
    return str(value)


def record(value, fields, label):
    """Validate one named record and serialize fields in the transport's declared order."""
    exact_keys(value,[k for k,_ in fields],label)
    return ' '.join(atom(value[k],kind) for k,kind in fields)


def rows(value):
    """Bound variable collection sizes before serializing or allocating model observations."""
    if not isinstance(value,list) or len(value)>1000000: raise ValueError('Invalid collection size')
    return value


def encode(document):
    """Validate native schema v1 and produce a typed in-memory stream for the reusable model.

    All observations/settings are supplied by this document. Estimated values
    must be strictly inside bounds and main deviations must already sum to zero.
    Parameter names and supported dimensions belong to this swordfish model,
    rather than a general arbitrary-species schema. No source files are opened.
    """
    exact_keys(document,['schema','units','model','parameters','recruitment_deviations'],'native input')
    if document['schema']!='quadra.swordfish.v1' or document['units']!=UNITS:
        raise ValueError('Unsupported native schema or units')
    m=document['model']
    exact_keys(m,['dimensions','biology','recruitment','blocks','population_length_bins','length_bins']+list(RECORDS)+list(TABLES)+[k for k,_ in SCALARS],'model')
    out=['QUADRA_SWORDFISH_NATIVE_V1',record(m['dimensions'],DIMENSIONS,'dimensions')]
    bio=m['biology'];exact_keys(bio,[k for k,_ in BIOLOGY]+['growth'],'biology')
    out.append(record({k:bio[k] for k,_ in BIOLOGY},BIOLOGY,'biology'))
    if len(rows(bio['growth']))!=2: raise ValueError('Exactly two sex-specific growth records required')
    out += [record(g,GROWTH,'growth') for g in bio['growth']]
    for name,fields in RECORDS.items():
        collection=rows(m[name]);out.append(str(len(collection)))
        out += [record(r,fields,name) for r in collection]
    for name in ['population_length_bins','length_bins']:
        collection=rows(m[name]);out.append(str(len(collection))+' '+' '.join(atom(v,'d') for v in collection))
    for name,(width,kind) in TABLES.items():
        collection=rows(m[name]);out.append(str(len(collection)))
        out += [atom(r,('a',width,kind)) for r in collection]
    blocks=rows(m['blocks']);out.append(str(len(blocks)))
    for block in blocks:
        block=rows(block)
        if len(block)%2: raise ValueError('Block years must be paired')
        out.append(str(len(block))+' '+' '.join(atom(v,'d') for v in block))
    out.append(record({k:m[k] for k,_ in SCALARS},SCALARS,'model scalars'))
    out.append(record(m['recruitment'],RECRUITMENT,'recruitment'))
    specs=rows(document['parameters']);out.append(str(len(specs)))
    names=set()
    for spec in specs:
        out.append(record(spec,PARAMETER,'parameter'))
        if spec['name'] in names: raise ValueError('Duplicate parameter')
        names.add(spec['name'])
        if spec['lower']>spec['upper'] or (spec['phase']>0 and not spec['lower']<spec['value']<spec['upper']):
            raise ValueError('Invalid parameter bounds/start: '+spec['name'])
        if spec['prior_type']!=0: raise ValueError('Only no-prior parameters supported')
        if not 0<=spec['block']<=len(blocks) or spec['block_function'] not in (0,1,2):
            raise ValueError('Unsupported parameter block')
    if names!=PARAMETER_NAMES: raise ValueError('Incomplete/unsupported parameter names')
    rec=m['recruitment']
    if not rec['early_start']<=rec['main_start']<=rec['main_end'] or rec['main_end']-rec['early_start']>1000:
        raise ValueError('Unsupported recruitment year range')
    dev=rows(document['recruitment_deviations']);out.append(str(len(dev)))
    years=set();main_sum=0.
    for r in dev:
        out.append(record(r,[('year','i'),('value','d')],'recruitment deviation'))
        if r['year'] in years: raise ValueError('Duplicate deviation year')
        years.add(r['year'])
        if not rec['minimum_deviation']<=r['value']<=rec['maximum_deviation']: raise ValueError('Deviation outside bounds')
        if r['year']>=rec['main_start']: main_sum+=r['value']
    if years!=set(range(rec['early_start'],rec['main_end']+1)) or abs(main_sum)>1e-8:
        raise ValueError('Incomplete deviation years or main deviations do not sum to zero')
    if not rec['early_start']<=rec['main_start']<=rec['main_end']: raise ValueError('Recruitment year ordering')
    if rec['cycles']!=0 or rec['read_deviations']!=0: raise ValueError('Unsupported recruitment options')
    validate_layout(document)
    return '\n'.join(out)+'\n'


def validate_layout(document):
    """Reject unsupported layouts and invalid observation references before launching C++."""
    m=document['model'];d=m['dimensions'];rec=m['recruitment']
    if (d['seasons'],d['subseasons'],d['plus_age'],d['areas'])!=(4,2,15,1):
        raise ValueError('Supported dimensions are four seasons, two subseasons, ages 0-15+, one area')
    if not 1<=d['first_year']<=d['last_year']<=9999 or d['last_year']-d['first_year']>500:
        raise ValueError('Unsupported model year range')
    if d['last_year']!=rec['main_end'] or rec['main_end']-rec['early_start']>1000:
        raise ValueError('Recruitment range must end at terminal year')
    if [f['id'] for f in m['fleets']]!=list(range(1,28)):
        raise ValueError('Fleet IDs must be ordered 1..27')
    if any(f['type']!=(1 if f['id']<=19 else 3) or f['area']!=1 or f['catch_multiplier']!=0 for f in m['fleets']):
        raise ValueError('Unsupported fleet structure')
    for name in ['index_specs','size_patterns','age_patterns','length_controls']:
        if len(m[name])!=27:raise ValueError(name+' requires 27 fleet settings')
    if any(r!=[0,0,0,0] for r in m['age_patterns']) or m['f_method']!=3:
        raise ValueError('Unsupported age selectivity or catch-conditioning method')
    for name in ['population_length_bins','length_bins']:
        bins=m[name]
        if len(bins)!=51 or any(a>=b for a,b in zip(bins,bins[1:])):
            raise ValueError(name+' requires 51 increasing bins')
    if [r[0] for r in m['q_setup']]!=list(range(20,28)):
        raise ValueError('Catchability setup requires index fleets 20..27')
    for r in m['catches']:
        if not 1<=r['fleet']<=19 or not 1<=r['season']<=4 or r['value']<0 or r['log_se']<0:
            raise ValueError('Invalid catch observation')
        if r['year']!=-999 and not d['first_year']<=r['year']<=d['last_year']:
            raise ValueError('Catch year outside model')
    for r in m['indices']:
        if not 20<=r['fleet']<=27 or (r['use'] and (r['value']<=0 or r['log_se']<=0 or not d['first_year']<=r['year']<=d['last_year'])):
            raise ValueError('Invalid abundance index')
    for r in m['compositions']:
        if not 1<=r['fleet']<=27 or r['sex']!=0 or r['partition']!=0 or any(v<0 for v in r['values']):
            raise ValueError('Unsupported/invalid length composition')
        if r['use'] and (r['sample_size']<=0 or sum(r['values'])<=0 or not d['first_year']<=r['year']<=d['last_year']):
            raise ValueError('Invalid included composition')


def unique_object(pairs):
    """Reject duplicate JSON object keys before values can be silently replaced."""
    result={}
    for key,value in pairs:
        if key in result: raise ValueError('Duplicate JSON field: '+key)
        result[key]=value
    return result


def load(text):
    """Load strict JSON, rejecting nonfinite constants and ambiguous duplicate keys."""
    def bad_constant(value):
        """Reject JSON extensions for NaN and infinity."""
        raise ValueError('Nonfinite JSON constant: '+value)
    value=json.loads(text,object_pairs_hook=unique_object,parse_constant=bad_constant)
    return value,encode(value)
