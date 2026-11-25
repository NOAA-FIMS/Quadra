"""Exercise SS3 conversion through the public command and compiled executable."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT/'build/pifsc_swordfish'
sys.path.insert(0,str(ROOT/'scripts'))
from native_input import load
import test_ss_par


class ConverterCLITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_ss_par.ParameterFileTests.setUpClass()
        cls.parameters=test_ss_par.ParameterFileTests.text

    def entrypoints(self):
        commands=[[str(ROOT/'swordfish')]]
        if NATIVE.is_file():commands.append([str(NATIVE)])
        return commands

    def sources(self,root):
        data=root/"observations ; 'quoted'.dat";control=root/'control settings.ctl'
        shutil.copyfile(ROOT/'data/raw/swo2023_v004.dat',data)
        shutil.copyfile(ROOT/'data/raw/swo2023_v007.ctl',control)
        return data,control

    def test_help_and_missing_options(self):
        for command in self.entrypoints():
            with self.subTest(command=command):
                top=subprocess.run(command+['--help'],capture_output=True,text=True)
                self.assertEqual(top.returncode,0,top.stderr);self.assertIn('--convert-ss3',top.stdout)
                help_result=subprocess.run(command+['--convert-ss3','--help'],capture_output=True,text=True)
                self.assertEqual(help_result.returncode,0,help_result.stderr)
                for flag in ['--data','--control','--parameters','--output']:self.assertIn(flag,help_result.stdout)
                missing=subprocess.run(command+['--convert-ss3'],capture_output=True,text=True)
                self.assertEqual(missing.returncode,2);self.assertIn('--data',missing.stderr)
                self.assertNotIn('Backend usage',missing.stderr)

    def test_conversion_and_optional_parameters(self):
        for command in self.entrypoints():
            for use_parameters in [False,True]:
                with self.subTest(command=command,parameters=use_parameters),tempfile.TemporaryDirectory(prefix='ss3 conversion ') as directory:
                    root=Path(directory);data,control=self.sources(root);output=root/'native model.json'
                    args=['--data',str(data),'--control',str(control),'--output',str(output)]
                    if use_parameters:
                        parameters=root/'initial vector.par';parameters.write_text(self.parameters)
                        args+=['--parameters',str(parameters)]
                    result=subprocess.run(command+['--convert-ss3']+args,cwd=root,capture_output=True,text=True)
                    self.assertEqual(result.returncode,0,result.stderr)
                    document,_=load(output.read_text());provenance=json.loads(output.with_suffix('.provenance.json').read_text())
                    self.assertEqual(len(document['parameters']),74)
                    self.assertFalse(provenance['report_sso_used']);self.assertFalse(provenance['covariance_used'])
                    values={r['name']:r['value'] for r in document['parameters']}
                    self.assertEqual(values['L_at_Amin_Mal_GP_1'],83.4 if use_parameters else 83.2)
                    repeated=subprocess.run(command+['--convert-ss3']+args,cwd=root,capture_output=True,text=True)
                    self.assertEqual(repeated.returncode,1);self.assertIn('Refusing to replace',repeated.stderr)
                    self.assertNotIn('Traceback',repeated.stderr)

    @unittest.skipUnless(NATIVE.is_file(),'Build with make pifsc_swordfish for compiled CLI checks')
    def test_relocated_binary_converts_without_python_or_scripts(self):
        with tempfile.TemporaryDirectory(prefix="standalone 'quoted' ") as directory:
            root=Path(directory);binary=root/'pifsc_swordfish'
            shutil.copy2(NATIVE,binary)
            data,control=self.sources(root);output=root/'native model.json'
            env=dict(os.environ,PATH='',SWORDFISH_PYTHON='/missing/python',SWORDFISH_CONVERTER='/missing/converter.py')
            result=subprocess.run([str(binary),'--convert-ss3','--data',str(data),'--control',str(control),'--output',str(output)],
                                  cwd=root,env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            load(output.read_text())
            self.assertEqual(sorted(p.name for p in root.iterdir()),sorted([binary.name,data.name,control.name,output.name,output.with_suffix('.provenance.json').name]))

    @unittest.skipUnless(NATIVE.is_file(),'Build with make pifsc_swordfish for compiled CLI checks')
    def test_native_found_through_path_from_another_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            env=dict(os.environ,PATH=str(NATIVE.parent)+os.pathsep+os.environ.get('PATH',''))
            result=subprocess.run([NATIVE.name,'--convert-ss3','--help'],cwd=directory,env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr);self.assertIn('--data',result.stdout)


if __name__=='__main__':unittest.main()
