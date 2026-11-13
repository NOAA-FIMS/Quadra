"""Run both control-only initialization experiments in a new scratch directory."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
from runtime import atomic_json, file_hash
from compare_control_start import compare
ROOT=Path(__file__).resolve().parents[1]


def main():
    """Execute fresh strict/matched-fixed fits, then read the saved reference only for comparison.

    Optimization receives no Report.sso, fitted parameter CSV or saved covariance.
    The matched-fixed sensitivity changes only the explicit fixed male-growth
    setting. Existing output directories are rejected rather than overwritten.
    Optimizer stdout/stderr is streamed to the terminal and flushed to a log.
    """
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='New directory; default build/control_start/cold-UTC')
    args=parser.parse_args()
    output=(args.output or ROOT/'build/control_start'/('cold-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))).resolve()
    output.mkdir(parents=True,exist_ok=False)
    identity={name:file_hash(ROOT/'data/raw'/name) for name in ['swo2023_v004.dat','swo2023_v007.ctl']}
    binary=ROOT/'build/swordfish_control_start'
    print('Run directory: '+str(output),flush=True)
    try:
        for name,mode in [('strict','strict'),('matched_fixed','matched-fixed')]:
            print('Starting control-only '+mode,flush=True)
            with (output/(name+'.log')).open('w') as stream:
                command=[str(binary),mode,str(output/name),'200']
                with subprocess.Popen(command,cwd=ROOT,stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT,text=True,bufsize=1) as process:
                    for line in process.stdout:
                        sys.stdout.write(line);sys.stdout.flush()
                        stream.write(line);stream.flush()
                    status=process.wait()
                    if status:
                        raise subprocess.CalledProcessError(status,command)
        if identity!={name:file_hash(ROOT/'data/raw'/name) for name in identity}:
            raise ValueError('Control inputs changed during fitting')
        compare(output/'strict',output/'matched_fixed',output/'comparison')
        atomic_json(output/'run_status.json',dict(status='passed',source_inputs_sha256=identity,binary_sha256=file_hash(binary),optimization_inputs='original dat/ctl factory; no fitted report/covariance'))
    except BaseException as error:
        atomic_json(output/'run_status.json',dict(status='failed',error=str(error),source_inputs_sha256=identity))
        raise
    print('Completed independent control-start experiments: '+str(output))


if __name__=='__main__':main()
