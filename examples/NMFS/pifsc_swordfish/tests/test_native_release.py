"""Native release safeguards and deterministic perturbations, without SS3 sources."""
import copy
import csv
from pathlib import Path
import sys
import tempfile
import subprocess
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from test_native_input import fixture
import run_native_release as release


class NativeReleaseTests(unittest.TestCase):
    """Check seed identity, covariance perturbations, source isolation and failed publication paths."""
    def test_seed_and_perturbations(self):
        """Preserve fixed settings, perturb reproducibly and reconstruct the dependent deviation."""
        model=fixture()
        next(p for p in model['parameters'] if p['name']=='SR_LN(R0)')['phase']=1
        expected={p['name']:p['value'] for p in model['parameters']}
        expected.update({f"RecruitmentDev_{d['year']}":d['value'] for d in model['recruitment_deviations']})
        release.verify_seed(model,expected)
        changed=copy.deepcopy(model);changed['parameters'][0]['value']+=1
        with self.assertRaises(ValueError):release.verify_seed(changed,expected)
        schema=release.coordinate_schema(model);factor=[[1. if i==j else 0. for j in range(len(schema))] for i in range(len(schema))]
        a=release.perturbed_input(model,factor,101,.1);b=release.perturbed_input(model,factor,101,.1)
        self.assertEqual(a,b);self.assertNotEqual(a,model);self.assertEqual(a['model'],model['model'])
        for p,q in zip(model['parameters'],a['parameters']):
            if p['phase']<=0:self.assertEqual(p,q)
        self.assertAlmostEqual(sum(d['value'] for d in a['recruitment_deviations'] if d['year']>=1985),0.,places=12)

    def test_metric_identity_and_positive_definiteness(self):
        """Read a known SPD matrix; reject wrong labels and indefinite directions."""
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'metric.csv';path.write_text('quantity,a,b\na,4,1\nb,1,2\n')
            factor=release.metric_factor(path,['a','b']);self.assertAlmostEqual(factor[0][0],2.)
            self.assertAlmostEqual(factor[1][0],.5)
            with self.assertRaises(ValueError):release.metric_factor(path,['b','a'])
            path.write_text('quantity,a,b\na,1,2\nb,2,1\n')
            with self.assertRaises(ValueError):release.metric_factor(path,['a','b'])

    def test_pdf_interpreter_keeps_virtual_environment(self):
        """Preserve the selected venv launcher rather than resolving it to the base Python."""
        with tempfile.TemporaryDirectory() as directory:
            interpreter=Path(directory)/'python';interpreter.symlink_to(sys.executable)
            with patch.object(sys,'argv',['release','--input','missing.json','--pdf-python',str(interpreter)]),patch.object(release,'release') as execute:
                release.main()
                self.assertEqual(execute.call_args.args[3],interpreter.absolute())

    def test_existing_destination_is_preserved(self):
        """Reject an existing release before reading missing inputs or starting child processes."""
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);out=root/'release';out.mkdir();(out/'keep').write_text('original')
            with self.assertRaises(FileExistsError):release.release(root,root/'missing.json',out,Path(sys.executable))
            self.assertEqual((out/'keep').read_text(),'original')

    def test_stage_failure_and_timeout(self):
        """Retain error logs/counters and terminate a timed-out process group."""
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);records=[]
            with self.assertRaises(RuntimeError):release.run_stage('failure',[sys.executable,'-c','print("recorded"); raise SystemExit(7)'],root,root,records)
            self.assertEqual(records[0]['returncode'],7)
            self.assertIn('recorded',(root/records[0]['log']).read_text())
            with self.assertRaises(subprocess.TimeoutExpired):release.run_stage('timeout',[sys.executable,'-c','import time; time.sleep(5)'],root,root,records,timeout=.05)
            self.assertNotEqual(records[-1]['returncode'],0)

    def test_snapshot_excludes_raw_inputs_and_old_outputs(self):
        """Copy source/targets into the correct layout without observation files or existing binaries."""
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);repo=root/'repo';case=repo/'examples/NMFS/pifsc_swordfish';case.mkdir(parents=True)
            for name in ['core','model','math','include','external/eigen','external/LBFGSpp','tests']:
                path=repo/name;path.mkdir(parents=True,exist_ok=True);(path/'source.hpp').write_text('source')
            (repo/'tests/test.cpp').write_text('source');(repo/'Makefile').write_text('source')
            for name in ['include','src','scripts','tests','data/raw','build']:
                path=case/name;path.mkdir(parents=True);(path/'old').write_text('not copied' if name in ['data/raw','build'] else 'source')
            for name in ['Makefile','swordfish']:(case/name).write_text('source')
            for name in release.TARGETS:
                path=case/'reference'/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('target')
            bundle=root/'bundle';frozen=release.copy_snapshot(case,bundle)
            self.assertFalse((frozen/'data/raw').exists());self.assertFalse((frozen/'build').exists())
            self.assertTrue((frozen/'src/old').is_file());self.assertTrue((bundle/'workspace/external/eigen/source.hpp').is_file())


if __name__=='__main__':unittest.main()
