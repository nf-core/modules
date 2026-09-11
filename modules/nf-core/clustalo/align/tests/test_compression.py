"""Exercise the CLUSTALO_ALIGN output shell expression without Nextflow.

Clustalo and pigz are stand-ins: these tests cover process completion and error
propagation, not alignment or pigz correctness. Run directly with Python 3.
Set CLUSTALO_MODULE to compare the same tests against an unpatched main.nf.
"""

import gzip
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


MODULE = Path(os.environ.get("CLUSTALO_MODULE", Path(__file__).parents[1] / "main.nf"))
ALIGNMENT = b">first\nAAAA-A\n>second\nAAATAA\n"

CLUSTALO = r'''
import os
from pathlib import Path
import sys

args = sys.argv[1:]
if os.environ.get("PRODUCER_STATUS"):
    sys.exit(int(os.environ["PRODUCER_STATUS"]))
output = args[args.index("-o") + 1]
with open(output, "wb") as handle:
    handle.write(b">first\nAAAA-A\n>second\nAAATAA\n")
print("alignment progress message")
'''

PIGZ = r'''
import gzip
import os
from pathlib import Path
import sys
import time

data = sys.stdin.buffer.read()
if os.environ.get("CONSUMER_STATUS"):
    sys.exit(int(os.environ["CONSUMER_STATUS"]))
if os.environ.get("DELAY_CONSUMER"):
    time.sleep(0.2)
sys.stdout.buffer.write(gzip.compress(data, mtime=0))
sys.stdout.buffer.flush()
Path("compression-complete").touch()
'''


class CompressionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        for name, source in (("clustalo", CLUSTALO), ("pigz", PIGZ)):
            executable = self.work / name
            executable.write_text("#!" + sys.executable + "\n" + source)
            executable.chmod(0o755)

    def run_output(self, compress=True, **environment):
        match = re.search(
            r'def write_output = compress \? "((?:\\.|[^"\\])*)" : "((?:\\.|[^"\\])*)"',
            MODULE.read_text(),
        )
        self.assertIsNotNone(match, "could not locate the module output expression")
        expression = match.group(1 if compress else 2)
        expression = expression.replace("${task.cpus}", "1").replace("${prefix}", "test")
        expression = expression.replace(r"\$", "$")
        env = dict(os.environ, PATH=str(self.work) + os.pathsep + os.environ["PATH"])
        for key in ("PRODUCER_STATUS", "CONSUMER_STATUS", "DELAY_CONSUMER"):
            env.pop(key, None)
        env.update(environment)
        # Files, not capture_output: inherited pipes can make subprocess wait for
        # an unwaited process substitution and accidentally hide the original bug.
        with (self.work / "stdout").open("w") as out, (self.work / "stderr").open("w") as err:
            return subprocess.run(
                ["bash", "-euo", "pipefail", "-c", "clustalo " + expression],
                cwd=self.work, env=env, stdout=out, stderr=err, timeout=10,
            ).returncode

    def test_compressed_output_is_complete(self):
        self.assertEqual(self.run_output(DELAY_CONSUMER="1"), 0)
        self.assertTrue((self.work / "compression-complete").exists())
        self.assertEqual(gzip.decompress((self.work / "test.aln.gz").read_bytes()), ALIGNMENT)

    def test_compressor_failure_is_propagated(self):
        # Drain all input before failing, so the producer cannot detect this via SIGPIPE.
        self.assertEqual(self.run_output(CONSUMER_STATUS="42"), 42)

    def test_producer_failure_is_preserved(self):
        self.assertEqual(self.run_output(PRODUCER_STATUS="17"), 17)

    def test_uncompressed_output_is_unchanged(self):
        self.assertEqual(self.run_output(compress=False, CONSUMER_STATUS="42"), 0)
        self.assertEqual((self.work / "test.aln").read_bytes(), ALIGNMENT)
        self.assertFalse((self.work / "test.aln.gz").exists())

    def test_stdout_messages_do_not_enter_alignment(self):
        self.assertEqual(self.run_output(), 0)
        self.assertIn("alignment progress message", (self.work / "stdout").read_text())
        self.assertEqual(gzip.decompress((self.work / "test.aln.gz").read_bytes()), ALIGNMENT)


if __name__ == "__main__":
    unittest.main()
