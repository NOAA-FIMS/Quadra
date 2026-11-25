"""Portable native-schema checks; these fixtures require no SS3 source files."""
import copy
import json
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from native_input import DIMENSIONS,GROWTH,BIOLOGY,RECRUITMENT,SCALARS,TABLES,PARAMETER_NAMES,UNITS,atom,encode,load


def field_default(kind):
    """Make typed scalar/array placeholders for a format-only synthetic fixture."""
    if isinstance(kind,tuple):return [field_default(kind[2]) for _ in range(kind[1])]
    return 0 if kind=='i' else 0.


def fixture():
    """Create a complete native schema example for validation tests, without real observations."""
    dimensions=dict(zip([k for k,_ in DIMENSIONS],[1975,2021,4,2,15,1,5.,[3.,3.,3.,3.]]))
    biology={k:field_default(kind) for k,kind in BIOLOGY}
    biology['growth']=[{k:1. for k,_ in GROWTH} for _ in range(2)]
    model=dict(dimensions=dimensions,biology=biology,catches=[],indices=[],compositions=[],blocks=[],population_length_bins=list(range(51)),length_bins=list(range(51)))
    model['fleets']=[dict(id=f,type=1 if f<=19 else 3,timing=1,area=1,catch_units=1,catch_multiplier=0,name=f'Fleet {f}') for f in range(1,28)]
    for name,(width,kind) in TABLES.items():model[name]=[[0]*width for _ in range(27 if name in ['index_specs','size_patterns','age_patterns','length_controls'] else 0)]
    model['q_setup']=[[f,1,0,0,0,1] for f in range(20,28)]
    model.update({k:field_default(kind) for k,kind in SCALARS});model['f_method']=3
    model['recruitment']={k:field_default(kind) for k,kind in RECRUITMENT}
    model['recruitment'].update(main_start=1985,main_end=2021,early_start=1965,phase=3,early_phase=4,minimum_deviation=-5.,maximum_deviation=5.)
    return dict(schema='quadra.swordfish.v1',units=UNITS.copy(),model=model,
                parameters=[dict(name=n,lower=-1000.,upper=1000.,value=0.,phase=-1,prior_type=0,block=0,block_function=0) for n in sorted(PARAMETER_NAMES)],
                recruitment_deviations=[dict(year=y,value=0.) for y in range(1965,2022)])


class NativeInputTests(unittest.TestCase):
    """Reject ambiguous JSON and malformed model states before the backend runs."""
    def test_valid_transport_and_roundtrip(self):
        """Check named-schema roundtrip and private transport version without source conversion."""
        model=fixture();document,payload=load(json.dumps(model))
        self.assertEqual(document,model);self.assertEqual(payload,encode(model))
        self.assertTrue(payload.startswith('QUADRA_SWORDFISH_NATIVE_V1\n'))

    def test_json_duplicates_and_nonfinite(self):
        """Reject duplicate object keys, NaN and infinity at the JSON boundary."""
        for text in ['{"schema":1,"schema":2}','{"x":NaN}','{"x":Infinity}']:
            with self.assertRaises(ValueError):load(text)

    def test_types_shapes_and_unknown_fields(self):
        """Reject booleans as numbers, fractional integers, bad arrays and misspelled keys."""
        for value,kind in [(True,'d'),(1.5,'i'),(float('inf'),'d'),([1,2],('a',3,'d'))]:
            with self.assertRaises(ValueError):atom(value,kind)
        model=fixture();model['unknown']=1
        with self.assertRaises(ValueError):encode(model)

    def test_bounds_names_and_deviation_constraints(self):
        """Reject incomplete/duplicate names, invalid starts and dependent-deviation inconsistencies."""
        for mode in ['missing','duplicate','bounds','sum','years','fleet','dimensions']:
            model=copy.deepcopy(fixture())
            if mode=='missing':model['parameters'].pop()
            if mode=='duplicate':model['parameters'][-1]=model['parameters'][0]
            if mode=='bounds':model['parameters'][0].update(phase=1,value=1000.)
            if mode=='sum':model['recruitment_deviations'][-1]['value']=1.
            if mode=='years':model['recruitment_deviations'].pop()
            if mode=='fleet':model['model']['fleets'][0]['id']=2
            if mode=='dimensions':model['model']['dimensions']['seasons']=2
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError):encode(model)


if __name__=='__main__':unittest.main()
