#!/usr/bin/env python3
"""Preserve supplied SWO inputs and extract selected SS3 3.30.20 report tables.

This is a version-specific intake tool, not a general SS3 parser or a model fit.
No external packages are required. Raw inputs remain byte-identical.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
FILES = [
    'swo2023_v004.dat', 'swo2023_v007.ctl', 'starter.ss', 'forecast.ss',
    'Report.sso',
    'ISC23_ANNEX11-Stock_Assessment_Report_for_North_Pacific_Ocean_Swordfish_FINAL.pdf',
    'North Pacific Swordfish 2023 Stock Assessment - SS3 Files & Reference Points Summary.docx',
    'North Pacific Swordfish 2023 Stock Assessment - SS3 Files & Reference Points Summary.pdf',
]


def read_lines(path):
    # Some exports contain literal escaped CRLF as well as actual line endings.
    """Read source text into lines using the supplied report encoding conventions.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.

    Returns
    -------
    Normalized source-line list.

    Notes
    -----
    Section/table matching is version specific. Raw files are preserved byte for byte and referenced by hashes. Extracted report tables are reference data, not an alternative observation source for the model.
    """
    return path.read_text(encoding='utf-8-sig').replace('\\r\\n', '\n').splitlines()


def section(lines, name):
    """Find an exact SS3 section label and return its bounded source lines; reject missing sections.

    Parameters
    ----------
    lines
        Source report lines after the documented encoding/newline normalization.
    name
        Exact section, output, field or stage identifier; matching is case sensitive.

    Returns
    -------
    Lines bounded by the requested report section and the next section marker.

    Notes
    -----
    Section/table matching is version specific. Raw files are preserved byte for byte and referenced by hashes. Extracted report tables are reference data, not an alternative observation source for the model.
    """
    start = next(i for i, s in enumerate(lines)
                 if re.match(r'^' + re.escape(name) + r' report:\d+', s))
    end = next((i for i in range(start + 1, len(lines))
                if re.match(r'^\S+ report:\d+', lines[i])), len(lines))
    return lines[start + 1:end]


def extract(lines, name, header_prefix, row_pattern):
    """Extract a checked header and pattern-matched rows from one version-specific report section.

    Parameters
    ----------
    lines
        Source report lines after the documented encoding/newline normalization.
    name
        Exact section, output, field or stage identifier; matching is case sensitive.
    header_prefix
        Expected header prefix within the selected version-specific report section.
    row_pattern
        Regular expression selecting supported data rows after the header.

    Returns
    -------
    Checked column header and supported extracted data rows.

    Notes
    -----
    Section/table matching is version specific. Raw files are preserved byte for byte and referenced by hashes. Extracted report tables are reference data, not an alternative observation source for the model.
    """
    block = section(lines, name)
    h = next(i for i, s in enumerate(block) if s.startswith(header_prefix))
    header = block[h].split()
    if len(set(header)) != len(header):
        raise ValueError(f'{name}: duplicate column names')
    rows = []
    for s in block[h + 1:]:
        if not s.strip():
            if rows:
                break
            continue
        if s.startswith('#'):
            if rows:
                break
            continue
        if not re.match(row_pattern, s):
            break
        row = s.split()
        if len(row) > len(header):
            raise ValueError(f'{name}: {len(row)} values for {len(header)} columns: {s[:100]}')
        if len(row) < len(header) and name not in ('DERIVED_QUANTITIES', 'PARAMETERS'):
            raise ValueError(f'{name}: truncated row: {s[:100]}')
        rows.append(row + [''] * (len(header) - len(row)))
    if not rows:
        raise ValueError(f'{name}: empty table')
    return header, rows


def write_csv(path, header, rows):
    """Write a header and extracted rows to a reference CSV while preserving column order.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.
    header
        Ordered output column names.
    rows
        Ordered report records or wrapped document-table rows, depending on the owning function.

    Notes
    -----
    Section/table matching is version specific. Raw files are preserved byte for byte and referenced by hashes. Extracted report tables are reference data, not an alternative observation source for the model.
    """
    with path.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def main():
    """Execute the import ss3 command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    Section/table matching is version specific. Raw files are preserved byte for byte and referenced by hashes. Extracted report tables are reference data, not an alternative observation source for the model.
    """
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True, help='Directory containing the eight supplied files')
    args = ap.parse_args()
    # Validate before copying. Never silently replace a different reference run.
    raw = ROOT / 'data/raw'
    for name in FILES:
        source = args.source / name
        if not source.is_file():
            raise FileNotFoundError(source)
        target = raw / name
        if target.exists() and target.read_bytes() != source.read_bytes():
            raise ValueError(f'Refusing to overwrite different raw input: {target}')
    raw.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name in FILES:
        source, target = args.source / name, raw / name
        if source.resolve() != target.resolve():
            shutil.copyfile(source, target)
        data = target.read_bytes()
        manifest.append({'file': name, 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
    reference = ROOT / 'reference'
    normalized = ROOT / 'data/normalized'
    reference.mkdir(exist_ok=True)
    normalized.mkdir(exist_ok=True)
    report = read_lines(raw / 'Report.sso')
    if not report[0].startswith('#V3.30.20.00;'):
        raise ValueError('Unsupported SS3 report version; inspect schema before importing')
    specs = [
        ('CATCH', 'Fleet Fleet_Name', r'^\d+\s', 'catch_fits.csv'),
        ('LIKELIHOOD', 'Component logL*Lambda', r'^\S+\s+[-+\d]', 'likelihood_components.csv'),
        ('PARAMETERS', 'Num Label Value', r'^\d+\s', 'parameters.csv'),
        ('DERIVED_QUANTITIES', 'Label Value', r'^\S+\s+[-+\d]', 'derived_quantities.csv'),
        ('TIME_SERIES', 'Area Yr Era Seas', r'^\d+\s', 'seasonal_time_series.csv'),
        ('Kobe_Plot', 'Yr ', r'^\d{4}\s', 'annual_stock_status.csv'),
        ('INDEX_2', 'Fleet Fleet_name', r'^\d+\s', 'index_fits.csv'),
        ('FIT_LEN_COMPS', 'Fleet Fleet_Name', r'^\d+\s', 'length_fit_summary.csv'),
    ]
    counts = {}
    tables = {}
    for name, prefix, pattern, filename in specs:
        header, rows = extract(report, name, prefix, pattern)
        write_csv(normalized / filename, header, rows)
        counts[filename] = len(rows)
        tables[name] = [dict(zip(header, row)) for row in rows]
    # Parse the data file's structural header, retaining original catch unit codes.
    dat = read_lines(raw / 'swo2023_v004.dat')
    active = [s.split('#', 1)[0].split() for s in dat if s.split('#', 1)[0].strip()]
    nfleets = int(active[9][0])
    fleet_rows = active[10:10 + nfleets]
    if len(fleet_rows) != nfleets or any(len(r) != 6 for r in fleet_rows):
        raise ValueError('Unexpected fleet schema')
    write_csv(normalized / 'fleets.csv',
              ['fleet', 'type', 'timing', 'area', 'catch_units', 'catch_multiplier', 'name'],
              [[i, *row] for i, row in enumerate(fleet_rows, 1)])
    counts['fleets.csv'] = nfleets
    dims = dict(zip(['start_year','end_year','seasons'], [int(active[i][0]) for i in range(3)]))
    dims.update(season_months=list(map(int, active[3])), subseasons=int(active[4][0]),
                spawn_month=float(active[5][0]), sexes=int(active[6][0]),
                plus_group_age=int(active[7][0]), areas=int(active[8][0]))
    if nfleets != 27:
        raise ValueError(f'Expected supplied 27-fleet assessment, got {nfleets}')
    historical = [r for r in tables['TIME_SERIES'] if r['Era'] == 'TIME']
    keys = [(r['Area'], r['Yr'], r['Seas']) for r in historical]
    expected = (dims['end_year'] - dims['start_year'] + 1) * dims['seasons'] * dims['areas']
    if len(keys) != expected or len(set(keys)) != expected:
        raise ValueError('Incomplete or duplicate historical seasonal rows')
    derived = {r['Label']: float(r['Value']) for r in tables['DERIVED_QUANTITIES']}
    # Independent cross-table consistency check, allowing printed precision.
    terminal = next(r for r in tables['Kobe_Plot'] if r['Yr'] == str(dims['end_year']))
    ratio = derived[f"SSB_{dims['end_year']}"] / derived['SSB_MSY']
    if abs(ratio - float(terminal['B/Bmsy'])) > 0.0001:
        raise ValueError('Terminal SSB ratio differs between report tables')
    summary = {'ss3_version': report[0], 'dimensions': dims, 'table_rows': counts,
               'historical_seasonal_rows': len(historical),
               'missing_for_original_initialization': ['ss.par'],
               'status': 'SS3 benchmark extraction only; no Quadra fit performed'}
    (reference / 'source_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (reference / 'intake.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
