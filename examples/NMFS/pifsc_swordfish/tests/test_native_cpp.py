"""Standalone C++ conversion parity, strict intake and runtime provenance."""
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT/'build/pifsc_swordfish'
sys.path.insert(0,str(ROOT/'scripts'))
from support.ss3_oracle import convert_files
from test_native_input import fixture

@unittest.skipUnless(BINARY.is_file(), 'Build with make pifsc_swordfish first')
class NativeCPPTests(unittest.TestCase):
    def invoke(self, *args, **kwargs):
        return subprocess.run([str(BINARY),*map(str,args)],capture_output=True,text=True,timeout=180,**kwargs)

    def test_invalid_cli_before_input_access(self):
        cases = [(['--input','missing','--output','unused','--mode','evaluate','--uncertainty'],'require --mode fit'),
                 (['--input','missing','--output','unused','--mode','evaluate','--optimizer-metric','missing'],'require --mode fit'),
                 (['--input','missing','--output','unused','--max-iterations','0'],'positive'),
                 (['--input','missing','--output','unused','--max-iterations','1.5'],'optimizer'),
                 (['--input','missing','--output','unused','--gradient-tolerance','NaN'],'optimizer'),
                 (['--input','missing','--output','unused','--input','other'],'Duplicate'),
                 (['--input','missing','--output','unused','--typo'],'Unknown'),
                 (['--convert-ss3','--data','a','--control','b','--input','c'],'Unknown')]
        for args, phrase in cases:
            with self.subTest(args=args):
                result=self.invoke(*args)
                self.assertEqual(result.returncode,2,result.stderr)
                self.assertIn(phrase.lower(),result.stderr.lower())

    def test_invalid_native_json_produces_no_output(self):
        mutations=[]
        for kind in ['unknown','type','names','bounds','sum','years','layout','array','units','block','rec_type']:
            doc=copy.deepcopy(fixture())
            if kind=='unknown':doc['model']['unknown']=0
            if kind=='type':doc['model']['dimensions']['first_year']=True
            if kind=='names':doc['parameters'].pop()
            if kind=='bounds':doc['parameters'][0].update(phase=1,value=1000.)
            if kind=='sum':doc['recruitment_deviations'][-1]['value']=1.
            if kind=='years':doc['recruitment_deviations'].pop()
            if kind=='layout':doc['model']['fleets'][0]['id']=2
            if kind=='array':doc['model']['biology']['natural_mortality'].pop()
            if kind=='units':doc['units']['weight']='lbs'
            if kind=='block':doc['parameters'][0]['block']=99
            if kind=='rec_type':doc['model']['recruitment']['early_start']=1965.5
            mutations.append(json.dumps(doc))
        mutations += ['{"schema":1,"schema":2}', '{"model":{"x":1,"x":2}}', '{"x":NaN}', '{"x":Infinity}', '{"x":1e999}', '{', '[]', '{"x":true} trailing']
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'bad.json';output=root/'must not exist'
            for text in mutations:
                with self.subTest(text=text[:50]):
                    source.write_text(text)
                    result=self.invoke('--input',source,'--output',output,'--mode','evaluate')
                    self.assertEqual(result.returncode,1,result.stderr)
                    self.assertFalse(output.exists())

    @unittest.skipUnless((ROOT/'data/raw/ss.par').is_file(),'Requires local SS3 fixtures')
    def test_conversion_parity_and_isolated_evaluation(self):
        data=ROOT/'data/raw/swo2023_v004.dat';control=ROOT/'data/raw/swo2023_v007.ctl';parameters=ROOT/'data/raw/ss.par'
        with tempfile.TemporaryDirectory(prefix="native C++ ' ") as directory:
            root=Path(directory);binary=root/'standalone';shutil.copy2(BINARY,binary)
            env=dict(os.environ,PATH='',SWORDFISH_PYTHON='/no/python',SWORDFISH_CONVERTER='/no/script')
            for use_parameters in [False,True]:
                source=root/f'cpp-{use_parameters}.json';expected=root/f'python-{use_parameters}.json'
                args=['--convert-ss3','--data',str(data),'--control',str(control),'--output',str(source)]
                if use_parameters:args+=['--parameters',str(parameters)]
                result=subprocess.run([str(binary),*args],env=env,cwd=root,capture_output=True,text=True,timeout=30)
                self.assertEqual(result.returncode,0,result.stderr)
                doc,side=convert_files(data,control,expected,parameters if use_parameters else None)
                self.assertEqual(json.loads(source.read_text()),doc)
                actual=json.loads(source.with_suffix('.provenance.json').read_text());reference=json.loads(side.read_text())
                self.assertEqual(actual.pop('native_json_sha256'),hashlib.sha256(source.read_bytes()).hexdigest());reference.pop('native_json_sha256')
                if use_parameters:
                    self.assertAlmostEqual(actual['saved_ss3_metadata'].pop('main_deviation_sum'),reference['saved_ss3_metadata'].pop('main_deviation_sum'),places=14)
                self.assertEqual(actual,reference)
                if use_parameters:
                    output=root/'evaluation output'
                    result=subprocess.run([str(binary),'--input',str(source),'--output',str(output),'--mode','evaluate'],env=env,cwd=root,capture_output=True,text=True,timeout=180)
                    self.assertEqual(result.returncode,0,result.stderr)
                    summary=json.loads((output/'quadra_fit_summary.json').read_text());status=json.loads((output/'run_status.json').read_text())
                    self.assertAlmostEqual(summary['objective'],289.299639305924,places=8)
                    self.assertEqual((output/'model_input.json').read_bytes(),source.read_bytes())
                    self.assertEqual(status['input_sha256'],hashlib.sha256(source.read_bytes()).hexdigest())
                    self.assertEqual(status['binary_sha256'],hashlib.sha256(binary.read_bytes()).hexdigest())
                    self.assertEqual(status['status'],'completed');self.assertFalse(status['ss3_files_used'])
                    before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()}
                    refused=self.invoke('--input',source,'--output',output,'--mode','evaluate')
                    self.assertNotEqual(refused.returncode,0)
                    self.assertEqual(before,{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()})

    @unittest.skipUnless((ROOT/'data/raw/ss.par').is_file(),'Requires local SS3 fixtures')
    def test_failed_fit_records_status(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'model.json';output=root/'incomplete'
            convert_files(ROOT/'data/raw/swo2023_v004.dat',ROOT/'data/raw/swo2023_v007.ctl',source,ROOT/'data/raw/ss.par')
            result=self.invoke('--input',source,'--output',output,'--max-iterations','1','--gradient-tolerance','1e-12')
            self.assertEqual(result.returncode,1,result.stderr)
            status=json.loads((output/'run_status.json').read_text())
            self.assertEqual(status['status'],'failed');self.assertEqual(status['returncode'],1)
            self.assertEqual((output/'model_input.json').read_bytes(),source.read_bytes())
            self.assertFalse((output/'quadra_fit_summary.json').exists())

    @unittest.skipUnless((ROOT/'data/raw/swo2023_v004.dat').is_file(),'Requires local SS3 fixtures')
    def test_ss3_unsupported_rows_truncation_and_duplicate_catches(self):
        raw=(ROOT/'data/raw/swo2023_v004.dat').read_text().replace('\\r\\n','\n')
        lines=raw.splitlines()
        catch_line=next(i for i,line in enumerate(lines) if len(line.partition('#')[0].split())==5)
        duplicate=lines[:catch_line]+[lines[catch_line]]+lines[catch_line:]
        texts=[raw.replace('0 #_N_age_bins','1 #_N_age_bins'),raw.rsplit('999',1)[0],
               raw.replace('10 15 20 25 30 35 40 45 50 55','11 15 20 25 30 35 40 45 50 55'),
               '\n'.join(duplicate),raw+'\n1\n']
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);data=root/'bad.dat';output=root/'model.json'
            for text in texts:
                with self.subTest(text=text[-50:]):
                    self.assertNotEqual(text,raw)
                    data.write_text(text)
                    result=self.invoke('--convert-ss3','--data',data,'--control',ROOT/'data/raw/swo2023_v007.ctl','--output',output)
                    self.assertEqual(result.returncode,1,result.stderr)
                    self.assertFalse(output.exists());self.assertFalse(output.with_suffix('.provenance.json').exists())

    @unittest.skipUnless((ROOT/'data/raw/ss.par').is_file(),'Requires local SS3 fixtures')
    def test_ss3_errors_and_existing_sidecar(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);output=root/'model.json';data=root/'bad.dat';control=ROOT/'data/raw/swo2023_v007.ctl'
            data.write_text((ROOT/'data/raw/swo2023_v004.dat').read_text().replace('#V3.30.20.00;','#V3.30.99.00;',1))
            result=self.invoke('--convert-ss3','--data',data,'--control',control,'--output',output)
            self.assertEqual(result.returncode,1);self.assertIn('version',result.stderr);self.assertFalse(output.exists())
            data.write_bytes((ROOT/'data/raw/swo2023_v004.dat').read_bytes())
            par=root/'bad.par';par.write_text((ROOT/'data/raw/ss.par').read_text()+'\n# unexpected:\n1\n')
            result=self.invoke('--convert-ss3','--data',data,'--control',control,'--parameters',par,'--output',output)
            self.assertEqual(result.returncode,1);self.assertFalse(output.exists())
            side=output.with_suffix('.provenance.json');side.write_text('keep')
            result=self.invoke('--convert-ss3','--data',data,'--control',control,'--output',output)
            self.assertEqual(result.returncode,1);self.assertIn('Refusing to replace',result.stderr)
            self.assertEqual(side.read_text(),'keep');self.assertFalse(output.exists())

if __name__=='__main__': unittest.main()
