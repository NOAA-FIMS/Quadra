"""Resource-free entry points for the compiled command help."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

BINARY = Path(__file__).resolve().parents[1] / "build" / "pifsc_swordfish"

@unittest.skipUnless(BINARY.is_file(), "Build with make pifsc_swordfish first")
class NativeHelpTests(unittest.TestCase):
    def test_help_aliases_work_without_resources_or_input(self):
        with tempfile.TemporaryDirectory() as directory:
            env = dict(os.environ, SWORDFISH_PYTHON="/missing/python",
                       SWORDFISH_CONVERTER="/missing/converter.py")
            outputs = []
            for arguments in ([], ["-h"], ["--help"]):
                result = subprocess.run([str(BINARY), *arguments], input="",
                                        capture_output=True, text=True, timeout=10,
                                        cwd=directory, env=env)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stderr, "")
                outputs.append(result.stdout)
            self.assertEqual(outputs[0], outputs[1])
            self.assertEqual(outputs[1], outputs[2])
            for phrase in ("--convert-ss3", "--parameters", "--input MODEL.json",
                           "--mode fit|evaluate", "typed stream", "EXAMPLES",
                           "--max-iterations", "entirely in C++"):
                self.assertIn(phrase, outputs[0])
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_invalid_arguments_point_to_help(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "unused-output"
            result = subprocess.run([str(BINARY), str(output)], input="",
                                    capture_output=True, text=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--help", result.stderr)
            self.assertFalse(output.exists())

if __name__ == "__main__":
    unittest.main()
