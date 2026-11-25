"""Independent SS3 conversion oracle for tests; no production CLI."""
from contextlib import ExitStack
import json
import hashlib
import os
from pathlib import Path
from .native_export import convert
from runtime import file_hash
from ss_par import parse as parse_parameters
from swordfish_inputs import parse_data, parse_control


def write_pair(output, document, provenance):
    """Create JSON and provenance exclusively, rolling back owned files on ordinary errors.

    Existing outputs and sidecars are never replaced. Both documents are fully
    serialized before creation. This is not a multi-file atomic transaction
    across an uncatchable process termination.
    """
    sidecar=output.with_suffix('.provenance.json')
    for path in [output,sidecar]:
        if path.exists() or path.is_symlink():raise FileExistsError('Refusing to replace '+str(path))
    native_text=json.dumps(document,indent=2,allow_nan=False)+'\n'
    provenance=dict(provenance,native_json_sha256=hashlib.sha256(native_text.encode('utf-8')).hexdigest())
    texts=[native_text,json.dumps(provenance,indent=2,allow_nan=False)+'\n']
    output.parent.mkdir(parents=True,exist_ok=True);owned=[]
    try:
        with ExitStack() as stack:
            streams=[]
            for path in [output,sidecar]:
                stream=stack.enter_context(path.open('x',encoding='utf-8'));owned.append(path);streams.append(stream)
            for stream,text in zip(streams,texts):
                stream.write(text);stream.flush();os.fsync(stream.fileno())
    except BaseException:
        for path in owned:path.unlink(missing_ok=True)
        raise
    return sidecar


def convert_files(data_path, control_path, output, parameter_path=None):
    """Read explicit SS3 paths once for intake and write a complete self-contained native input.

    No hardcoded source filenames, Report.sso, forecast files, covariance or C++
    compiler are needed. Optional .par values override all native starts and
    historical deviations, including fixed settings. Otherwise the original
    control starts and zero deviations are used, without silent feasibility
    adjustments. Input hashes and saved .par metadata go into a separate sidecar.
    """
    data_path,control_path,output=Path(data_path),Path(control_path),Path(output)
    parameter_path=Path(parameter_path) if parameter_path is not None else None
    for path in [output,output.with_suffix('.provenance.json')]:
        if path.exists() or path.is_symlink():raise FileExistsError('Refusing to replace '+str(path))
    sources={'data':data_path,'control':control_path}
    if parameter_path is not None:sources['parameters']=parameter_path
    identities={role:dict(file=path.name,sha256=file_hash(path)) for role,path in sources.items()}
    observations=parse_data(data_path);controls=parse_control(control_path,observations)
    inputs=dict(data=observations,control=controls)
    values,metadata=parse_parameters(parameter_path,inputs) if parameter_path is not None else (None,None)
    document=convert(inputs,values)
    if any(file_hash(path)!=identities[role]['sha256'] for role,path in sources.items()):
        raise ValueError('SS3 inputs changed during conversion')
    provenance=dict(schema='quadra.ss3_to_pifsc_swordfish.provenance.v1',
                    native_schema=document['schema'],supported_source='NP swordfish SS3 3.30.20 assessment layout',
                    sources=identities,initialization='supplied ss.par values' if values is not None else 'control initial values; historical deviations zero',
                    saved_ss3_metadata=metadata,report_sso_used=False,covariance_used=False,
                    counts=dict(parameters=len(document['parameters']),historical_deviations=len(document['recruitment_deviations']),
                                catches=len(observations['catches']),indices=len(observations['indices']),
                                length_compositions=len(observations['length_compositions'])))
    sidecar=write_pair(output,document,provenance)
    return document,sidecar
