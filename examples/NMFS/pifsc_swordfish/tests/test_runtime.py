"""Failure-path checks for content-bound covariance and non-overwriting publication."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import runtime
import run_release

class RuntimeTests(unittest.TestCase):
    """Ensure corrupt or stale results cannot be published as successful releases."""
    def test_atomic_json_preserves_previous_value_on_failure(self):
        """Nonfinite JSON must leave the old artifact intact and remove its temporary file."""
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'status.json'
            runtime.atomic_json(path,{'value':1})
            with self.assertRaises(ValueError):runtime.atomic_json(path,{'value':float('nan')})
            self.assertEqual(json.loads(path.read_text()),{'value':1})
            self.assertEqual(len(list(Path(directory).iterdir())),1)

    def test_covariance_rejects_tampering_and_changed_sources(self):
        """Every stamped artifact and source identity must still match before reuse."""
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'build').mkdir()
            identity={'sources':{'model.cpp':'original'},'fitted_parameters':'fit'}
            with patch('runtime.covariance_inputs',return_value=identity):
                runtime.begin_covariance(root)
                for name in ['coordinate_covariance.csv','uncertainty_coordinates.csv','hessian_coarse.csv','hessian_fine.csv','curvature_summary.json']:(root/'build'/name).write_text('original')
                runtime.finish_covariance(root);runtime.verify_covariance(root)
                (root/'build/coordinate_covariance.csv').write_text('altered')
                with self.assertRaises(ValueError):runtime.verify_covariance(root)
                (root/'build/coordinate_covariance.csv').write_text('original')
                identity['sources']['model.cpp']='changed'
                with self.assertRaises(ValueError):runtime.verify_covariance(root)

    def test_no_existing_destination_is_replaced(self):
        """An existing destination is rejected before lock creation or source copying."""
        with tempfile.TemporaryDirectory() as directory:
            destination=Path(directory)/'release';destination.mkdir();(destination/'keep').write_text('original')
            with self.assertRaises(FileExistsError):run_release.release(Path(directory),destination)
            self.assertEqual((destination/'keep').read_text(),'original')

    def test_failed_run_is_retained_but_not_published(self):
        """A preflight failure releases the lock and records a failed, unpublished bundle."""
        with tempfile.TemporaryDirectory() as directory:
            destination=Path(directory)/'release'
            with patch('run_release.verify_sources',side_effect=ValueError('bad source')):
                with self.assertRaises(ValueError):run_release.release(Path(directory),destination)
            self.assertFalse(destination.exists());self.assertFalse((Path(directory)/'.release.lock').exists())
            failures=list(Path(directory).glob('.release.failed-*'))
            self.assertEqual(len(failures),1)
            self.assertEqual(json.loads((failures[0]/'failure.json').read_text())['status'],'failed')

    def test_stage_timeout_records_failure(self):
        """A timed-out stage must stop its process group and retain a failure log record."""
        import subprocess
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);records=[]
            with self.assertRaises(subprocess.TimeoutExpired):
                run_release.run_stage('timeout',[sys.executable,'-c','import time; time.sleep(60)'],root,root,records,timeout=.1)
            self.assertEqual(len(records),1)
            self.assertTrue(records[0]['interrupted'])
            self.assertNotEqual(records[0]['returncode'],0)
