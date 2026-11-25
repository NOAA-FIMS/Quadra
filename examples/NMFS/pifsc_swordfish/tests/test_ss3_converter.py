"""Check the supported SS3 intake utility without requiring Jon's original .par file."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from support.ss3_oracle import convert_files, write_pair
from native_input import load
import test_ss_par


class ConverterTests(unittest.TestCase):
    """Verify explicit filenames, initialization provenance and non-overwriting output behavior."""
    @classmethod
    def setUpClass(cls):
        """Reuse the synthetic parameter-file fixture; no fitted report or downloaded .par is read."""
        test_ss_par.ParameterFileTests.setUpClass()
        cls.parameter_text=test_ss_par.ParameterFileTests.text

    def sources(self, directory):
        """Copy only observation/control files under arbitrary names, without Report.sso."""
        directory=Path(directory);data=directory/'observations.input';control=directory/'settings.input'
        shutil.copyfile(ROOT/'data/raw/swo2023_v004.dat',data)
        shutil.copyfile(ROOT/'data/raw/swo2023_v007.ctl',control)
        return data,control

    def test_control_initialization(self):
        """Preserve control starts/zero deviations and bind native output to source hashes."""
        with tempfile.TemporaryDirectory() as directory:
            data,control=self.sources(directory);output=Path(directory)/'model.json'
            document,sidecar=convert_files(data,control,output)
            parsed,_=load(output.read_text());self.assertEqual(document,parsed)
            values={r['name']:r['value'] for r in document['parameters']}
            self.assertEqual(values['SR_LN(R0)'],4.5);self.assertEqual(values['L_at_Amin_Mal_GP_1'],83.2)
            self.assertTrue(all(r['value']==0 for r in document['recruitment_deviations']))
            provenance=json.loads(sidecar.read_text())
            self.assertFalse(provenance['report_sso_used']);self.assertIsNone(provenance['saved_ss3_metadata'])
            self.assertEqual(provenance['sources']['data']['file'],'observations.input')
            self.assertEqual(provenance['sources']['data']['sha256'],hashlib.sha256(data.read_bytes()).hexdigest())
            self.assertEqual(provenance['native_json_sha256'],hashlib.sha256(output.read_bytes()).hexdigest())

    def test_parameter_initialization(self):
        """Import fixed settings and deviations directly from a supplied parameter snapshot."""
        with tempfile.TemporaryDirectory() as directory:
            data,control=self.sources(directory);parameters=Path(directory)/'initial.values'
            parameters.write_text(self.parameter_text);output=Path(directory)/'model.json'
            document,sidecar=convert_files(data,control,output,parameters)
            values={r['name']:r['value'] for r in document['parameters']}
            self.assertEqual(values['L_at_Amin_Mal_GP_1'],83.4)
            self.assertEqual(document['recruitment_deviations'][0]['value'],.01)
            self.assertEqual(len(document['parameters']),74);self.assertEqual(len(document['recruitment_deviations']),57)
            provenance=json.loads(sidecar.read_text());self.assertEqual(provenance['saved_ss3_metadata']['ss3_parameter_count'],86)

    def test_existing_outputs_and_invalid_inputs(self):
        """Protect either existing output, and leave no files after a malformed input fails."""
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);data,control=self.sources(root);output=root/'model.json';sidecar=output.with_suffix('.provenance.json')
            sidecar.write_text('keep')
            with self.assertRaises(FileExistsError):convert_files(data,control,output)
            self.assertEqual(sidecar.read_text(),'keep');self.assertFalse(output.exists());sidecar.unlink()
            output.write_text('keep')
            with self.assertRaises(FileExistsError):convert_files(data,control,output)
            self.assertEqual(output.read_text(),'keep');output.unlink()
            data.write_text('malformed')
            with self.assertRaises(ValueError):convert_files(data,control,output)
            self.assertFalse(output.exists());self.assertFalse(sidecar.exists())

    def test_pair_rollback(self):
        """A sidecar creation failure removes only the newly created native output."""
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'model.json';sidecar=output.with_suffix('.provenance.json');original=Path.open
            def fail_sidecar(path,*args,**kwargs):
                """Simulate an I/O failure after reserving the first output."""
                if path==sidecar:raise OSError('sidecar unavailable')
                return original(path,*args,**kwargs)
            with patch.object(Path,'open',new=fail_sidecar):
                with self.assertRaises(OSError):write_pair(output,{},dict(source='test'))
            self.assertFalse(output.exists());self.assertFalse(sidecar.exists())



if __name__=='__main__':unittest.main()
