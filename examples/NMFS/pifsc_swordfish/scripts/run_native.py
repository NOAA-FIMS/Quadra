"""Run normalized JSON through the reusable native Quadra model; no SS3 dependencies."""
import argparse
import hashlib
from pathlib import Path
import subprocess
from native_input import load
from runtime import atomic_json, file_hash

ROOT=Path(__file__).resolve().parents[1]


def main():
    """Validate JSON, stream native data to the backend and save exact input provenance.

    The binary is compiled once by make native-swordfish. Dataset changes need
    no compilation, generated C++ factory, SS3 source files or reference outputs.
    Existing output directories are rejected. Progress streams to the terminal.
    """
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--mode',choices=['fit','evaluate'],default='fit')
    parser.add_argument('--max-iterations',type=int,default=200)
    parser.add_argument('--gradient-tolerance',type=float,default=1e-4)
    parser.add_argument('--uncertainty',action='store_true',help='Fresh local covariance and reference-point SEs/CVs/95%% intervals (fit mode)')
    parser.add_argument('--optimizer-metric',type=Path,help='Optional labeled positive-definite covariance used only for optimizer scaling')
    args=parser.parse_args();output=args.output.resolve()
    if args.optimizer_metric and args.mode!='fit':parser.error('--optimizer-metric requires --mode fit')
    if args.uncertainty and args.mode!='fit':parser.error('--uncertainty requires --mode fit')
    if output.exists():raise FileExistsError('Output directory already exists: '+str(output))
    binary=ROOT/'build/swordfish_native'
    if not binary.is_file():raise FileNotFoundError('Build the native runner first: make native-swordfish')
    source=args.input.read_bytes();_,payload=load(source.decode('utf-8'))
    print('Native input: '+str(args.input.resolve()),flush=True)
    command=[str(binary),str(output),args.mode,str(args.max_iterations),str(args.gradient_tolerance)]
    if args.uncertainty:command.append('uncertainty')
    metric_hash=None
    if args.optimizer_metric:
        metric_hash=file_hash(args.optimizer_metric)
        command.append('metric='+str(args.optimizer_metric.resolve()))
    result=subprocess.run(command,input=payload,text=True)
    metric_unchanged=not args.optimizer_metric or file_hash(args.optimizer_metric)==metric_hash
    if output.is_dir():
        (output/'model_input.json').write_bytes(source)
        atomic_json(output/'run_status.json',dict(status='completed' if result.returncode==0 and metric_unchanged else 'failed',mode=args.mode,
                    input_sha256=hashlib.sha256(source).hexdigest(),binary_sha256=file_hash(binary),ss3_files_used=False,uncertainty_requested=args.uncertainty,optimizer_metric_sha256=metric_hash,optimizer_metric_unchanged=metric_unchanged,returncode=result.returncode))
    result.check_returncode()
    if not metric_unchanged:raise ValueError('Optimizer metric changed during run')


if __name__=='__main__':main()
