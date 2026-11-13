"""Read the supported NP swordfish SS3 .par layout without using Report.sso."""
import math
import re


def parse(path, inputs):
    """Map SS3 parameter groups in control-file order to all 131 native values.

    The mapping covers 20 biology, 5 stock-recruit, 1 initial-F, 8 Q and
    40 selectivity values (including expanded time blocks), plus 20 early and
    37 main recruitment deviations. Dummy, checksum and forecast entries are
    recorded separately. Unknown, duplicate, missing, nonfinite or incorrectly
    sized groups are rejected. The .par header describes its saved SS3 fit;
    it does not prove that these were the original pre-optimization values.
    """
    text = path.read_text()
    if not text.strip():
        raise ValueError('Empty SS3 parameter file')
    header = re.fullmatch(
        r'# Number of parameters = (\d+) Objective function value = (\S+)\s+'
        r'Maximum gradient component = (\S+)\s*', text.splitlines()[0])
    if not header:
        raise ValueError('Unsupported SS3 parameter header')
    metadata = dict(ss3_parameter_count=int(header[1]), ss3_objective=float(header[2]),
                    ss3_maximum_gradient=float(header[3]))
    if metadata['ss3_parameter_count'] != 86 or not all(math.isfinite(metadata[k]) for k in ['ss3_objective', 'ss3_maximum_gradient']):
        raise ValueError('Unsupported parameter count or nonfinite SS3 metadata')
    groups = {}
    current = None
    for line in text.splitlines()[1:]:
        if not line.strip():
            continue
        if line.startswith('# '):
            match = re.fullmatch(r'# ([^:]+):\s*', line)
            if not match or match[1] in groups:
                raise ValueError('Invalid or duplicate SS3 parameter group')
            current = match[1]
            groups[current] = []
        else:
            if current is None:
                raise ValueError('Values outside a parameter group')
            values = [float(v) for v in line.split()]
            if not all(math.isfinite(v) for v in values):
                raise ValueError('Nonfinite SS3 parameter')
            groups[current].extend(values)
    control = inputs['control']
    specs = control['parameters']
    families = [('MGparm', specs[:20]), ('SR_parm', specs[20:25]),
                ('init_F', [p for p in specs if p['name'].startswith('InitF_')]),
                ('Q_parm', [p for p in specs if p['name'].startswith('LnQ_')]),
                ('selparm', [p for p in specs if p['name'].startswith('SizeSel_')])]
    if [len(rows) for _, rows in families] != [20, 5, 1, 8, 40]:
        raise ValueError('Unsupported control parameter layout')
    expected = {'recdev_early': 20, 'recdev1': 37, 'Fcast_recruitments': 1,
                'dummy_parm': 1, 'checksum999': 1}
    native = {}
    for family, rows in families:
        for i, spec in enumerate(rows, 1):
            key = f'{family}[{i}]'
            expected[key] = 1
            if key not in groups or len(groups[key]) != 1:
                raise ValueError('Missing or incorrectly sized group: ' + key)
            value = groups[key][0]
            lo, hi = spec['columns'][:2]
            if spec['columns'][6] > 0 and not lo < value < hi:
                raise ValueError('Estimated parameter outside open bounds: ' + spec['name'])
            native[spec['name']] = value
    if groups.keys() != expected.keys() or any(len(groups[k]) != n for k, n in expected.items()):
        raise ValueError('Unexpected, missing or incorrectly sized SS3 groups')
    if groups['dummy_parm'] != [1.] or groups['Fcast_recruitments'] != [0.]:
        raise ValueError('Unsupported dummy or forecast recruitment setting')
    rec = control['recruitment']
    if (int(rec['advanced'][0]), rec['start'], rec['end']) != (1965, 1985, 2021):
        raise ValueError('Unsupported recruitment year mapping')
    if abs(sum(groups['recdev1'])) > 1e-8:
        raise ValueError('Main recruitment deviations do not sum to zero')
    for first, key in [(1965, 'recdev_early'), (1985, 'recdev1')]:
        for year, value in enumerate(groups[key], first):
            if not -5 <= value <= 5:
                raise ValueError('Recruitment deviation outside model bounds')
            native[f'RecruitmentDev_{year}'] = value
    metadata.update(checksum999=groups['checksum999'][0],
                    forecast_recruitment=groups['Fcast_recruitments'][0],
                    dummy_parameter=groups['dummy_parm'][0],
                    main_deviation_sum=sum(groups['recdev1']),
                    mapped_native_parameters=len(specs), historical_deviations=57)
    return native, metadata
