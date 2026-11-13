import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from swordfish_inputs import load_inputs, parse_data, validate_data


class InputTests(unittest.TestCase):
    """Verify the original input schema, exclusions and unsupported-layout rejection."""
    @classmethod
    def setUpClass(cls):
        """Verify setUpClass. This helper belongs to the source validation test suite."""
        cls.model=load_inputs(ROOT/'data/raw')

    def test_observations_and_exclusions(self):
        """Verify observations and exclusions. This helper belongs to the source validation test suite."""
        d=self.model['data']
        self.assertEqual(len(d['catches']),2223)
        self.assertEqual(len(d['indices']),164)
        self.assertEqual(len(d['length_compositions']),446)
        excluded=[r for r in d['length_compositions'] if not r['use']]
        self.assertEqual(len(excluded),95)
        self.assertTrue(all(r['year']<0 for r in excluded))
        self.assertTrue(all(len(r['values'])==102 for r in d['length_compositions']))
        # Sex=0 samples carry two input blocks; do not reinterpret these as
        # two independent samples or silently add both to sample size.
        self.assertEqual(d['length_compositions'][0]['sample_size'],17.6)
        self.assertEqual(d['length_compositions'][0]['values'][:51],
                         d['length_compositions'][0]['values'][51:])

    def test_special_index_and_weights(self):
        """Verify special index and weights. This helper belongs to the source validation test suite."""
        ctl=self.model['control']
        self.assertEqual(self.model['data']['index_specs'][24],[25,33,0,0])
        self.assertIn([1.,24.,0.11],ctl['variance_adjustments'])
        self.assertIn([1.,24.,1.,0.,0.],ctl['lambda_changes'])
        self.assertIn([1.,27.,1.,0.,0.],ctl['lambda_changes'])
        self.assertEqual(len(ctl['parameters']),74)
        p={x['name']:x['columns'][2] for x in ctl['parameters']}
        self.assertEqual(p['L_at_Amin_Mal_GP_1'],83.2)
        self.assertEqual(ctl['settlement'],[1.,7.,1.,0.])
        self.assertEqual(ctl['size_patterns'][4],[15,0,0,4])
        self.assertEqual(ctl['blocks'],[[2004.,2009.],[1994.,1997.]])

    def test_duplicate_catch_rejected(self):
        """Verify duplicate catch rejected. This helper belongs to the source validation test suite."""
        data=copy.deepcopy(self.model['data'])
        data['catches'].append(data['catches'][0])
        with self.assertRaisesRegex(ValueError,'Duplicate catch'):
            validate_data(data)

    def test_unsupported_data_and_truncation_rejected(self):
        """Verify unsupported data and truncation rejected. This helper belongs to the source validation test suite."""
        raw=(ROOT/'data/raw/swo2023_v004.dat').read_text().replace('\\r\\n','\n')
        for text in [raw.replace('0 #_N_age_bins','1 #_N_age_bins'), raw.rsplit('999',1)[0],
                     raw.replace('10 15 20 25 30 35 40 45 50 55','11 15 20 25 30 35 40 45 50 55')]:
            with tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/'invalid.dat'; path.write_text(text)
                with self.assertRaises(ValueError): parse_data(path)

    def test_no_nonfinite_observations(self):
        """Verify no nonfinite observations. This helper belongs to the source validation test suite."""
        data=copy.deepcopy(self.model['data'])
        data['indices'][0]['value']=float('nan')
        with self.assertRaises(ValueError): validate_data(data)


if __name__=='__main__': unittest.main()
