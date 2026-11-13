"""Check .par structural mapping without distributing the supplied parameter file."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from ss_par import parse
from swordfish_inputs import load_inputs


class ParameterFileTests(unittest.TestCase):
    """Exercise native group ordering, recruitment constraints and malformed input rejection."""
    @classmethod
    def setUpClass(cls):
        """Construct a synthetic .par from control starts; no Jon-provided file is needed."""
        cls.inputs = load_inputs(ROOT/'data/raw')
        specs = cls.inputs['control']['parameters']
        families = [('MGparm', specs[:20]), ('SR_parm', specs[20:25]),
                    ('init_F', [p for p in specs if p['name'].startswith('InitF_')]),
                    ('Q_parm', [p for p in specs if p['name'].startswith('LnQ_')]),
                    ('selparm', [p for p in specs if p['name'].startswith('SizeSel_')])]
        lines = ['# Number of parameters = 86 Objective function value = 289.3 Maximum gradient component = 0.00005', '# dummy_parm:', '1']
        for family, rows in families:
            for index, spec in enumerate(rows, 1):
                value = 83.4 if spec['name'] == 'L_at_Amin_Mal_GP_1' else spec['columns'][2]
                lines += [f'# {family}[{index}]:', str(value)]
        lines += ['# recdev_early:', ' '.join(['0.01']*20), '# recdev1:',
                  ' '.join(['0.1', '-0.1']+['0']*35), '# Fcast_recruitments:', '0', '# checksum999:', '6.9']
        cls.text = '\n'.join(lines)+'\n'

    def read(self, text):
        """Parse temporary text through the public file reader and actual control schema."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'ss.par'; path.write_text(text)
            return parse(path, self.inputs)

    def test_complete_mapping(self):
        """Check every native mapping, expanded selectivity block and deviation year."""
        native, metadata = self.read(self.text)
        self.assertEqual(len(native), 131)
        for spec in self.inputs['control']['parameters']:
            expected = 83.4 if spec['name'] == 'L_at_Amin_Mal_GP_1' else spec['columns'][2]
            self.assertEqual(native[spec['name']], expected)
        self.assertEqual(native['RecruitmentDev_1965'], .01)
        self.assertEqual(native['RecruitmentDev_1984'], .01)
        self.assertEqual(native['RecruitmentDev_1985'], .1)
        self.assertEqual(native['RecruitmentDev_1986'], -.1)
        self.assertEqual(native['RecruitmentDev_2021'], 0)
        self.assertEqual(metadata['ss3_parameter_count'], 86)
        self.assertEqual(metadata['main_deviation_sum'], 0)

    def test_malformed_groups(self):
        """Reject missing/duplicate/extra/nonfinite groups and nonzero recruitment sums."""
        for text in ['', self.text+'# dummy_parm:\n1\n', self.text+'# unknown:\n0\n',
                     self.text.replace('# selparm[40]:', '# selparm[41]:'),
                     self.text.replace('# dummy_parm:\n1', '# dummy_parm:\nnan'),
                     self.text.replace('# dummy_parm:\n1', '# dummy_parm:\n1 2'),
                     self.text.replace('0.1 -0.1', '0.1 0.1')]:
            with self.subTest(text=text[-100:]):
                with self.assertRaises(ValueError): self.read(text)


if __name__ == '__main__':
    unittest.main()
