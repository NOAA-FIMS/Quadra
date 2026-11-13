"""Strict, row-oriented reader for the supplied SS3 3.30.20 swordfish inputs.

Comments are provenance, never authoritative fleet mappings or instructions.
Unsupported layouts fail rather than silently dropping observations/features.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import math


@dataclass
class Row:
    """One active SS3 row with source provenance.

    Attributes
    ----------
    values : list
        Raw active tokens before numeric conversion. Comments are not tokens.
    line : int
        One-based physical source line used in validation errors.
    comment : str
        Text after the source comment marker, retained without interpreting it
        as a model setting or instruction.
    """
    values: list
    line: int
    comment: str


class Reader:
    """Stateful row-oriented reader for the supported SS3 version.

    Attributes
    ----------
    path : pathlib.Path
        Source file used to qualify every location-sensitive error.
    rows : list[Row]
        Active rows in source order, retaining original line numbers/comments.
    i : int
        Index of the next unconsumed active row. take/numbers/scalar/table
        advance this cursor; callers must not rewind it to ignore bad rows.

    Notes
    -----
    finish consumes the required 999 terminator and rejects remaining active
    rows. Parsing a prefix is not accepted as a successful input import.
    """
    def __init__(self, path):
        """Read a UTF-8 SS3 file, verify its version prefix and retain active token rows with original source locations.

        Parameters
        ----------
        path
            Filesystem input/output path in the schema described by the function summary.
        """
        self.path = Path(path)
        text = self.path.read_text(encoding='utf-8-sig').replace('\\r\\n', '\n')
        if not text.startswith('#V3.30.20.00;'):
            raise ValueError(f'{path}: unsupported version')
        self.rows = []
        for line, s in enumerate(text.splitlines(), 1):
            active, _, comment = s.partition('#')
            if active.strip():
                self.rows.append(Row(active.split(), line, comment.strip()))
        self.i = 0

    def take(self, n=None):
        """Consume one active row and optionally require an exact token count; fail with a source location on EOF or mismatch.

        Parameters
        ----------
        n
            Required row width, or None when the source function intentionally leaves width unconstrained.

        Returns
        -------
        Next active Row; the reader position advances exactly once.

        Notes
        -----
        Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
        """
        if self.i >= len(self.rows):
            raise ValueError(f'{self.path}: unexpected end of input')
        row = self.rows[self.i]
        self.i += 1
        if n is not None and len(row.values) != n:
            raise ValueError(f'{self.path}:{row.line}: expected {n} columns, got {len(row.values)}')
        return row

    def numbers(self, n=None):
        """Consume one numeric row and reject malformed or nonfinite values with the original file/line context.

        Parameters
        ----------
        n
            Required row width, or None when the source function intentionally leaves width unconstrained.

        Returns
        -------
        Finite float values from the consumed row.

        Notes
        -----
        Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
        """
        row = self.take(n)
        try:
            values = list(map(float, row.values))
            if not all(math.isfinite(v) for v in values):
                raise ValueError('nonfinite value')
            return values
        except ValueError as exc:
            raise ValueError(f'{self.path}:{row.line}: invalid numeric row') from exc

    def scalar(self):
        """Consume a row containing exactly one finite numeric value.

        Returns
        -------
        One finite numeric scalar.

        Notes
        -----
        Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
        """
        return self.numbers(1)[0]

    def integer(self):
        """Consume one scalar and require it to represent an exact integer.

        Returns
        -------
        One exactly integral source value.

        Notes
        -----
        Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
        """
        value = self.scalar()
        if value != int(value):
            raise ValueError(f'{self.path}: expected integer, got {value}')
        return int(value)

    def expect(self, values):
        """Require the next row to equal a supported setting exactly instead of silently accepting another SS3 configuration.

        Parameters
        ----------
        values
            Numeric source values whose required integer/setting constraints are checked.

        Notes
        -----
        Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
        """
        actual = self.numbers(len(values))
        if actual != values:
            raise ValueError(f'{self.path}: unsupported setting {actual}, expected {values}')

    def table(self, n):
        """Read rows of a fixed-width source table until the supported sentinel is encountered.

        Parameters
        ----------
        n
            Required row width, or None when the source function intentionally leaves width unconstrained.

        Returns
        -------
        Numeric rows up to, but excluding, the supported table sentinel.
        """
        rows = []
        while True:
            row = self.numbers(n)
            if row[0] == -9999:
                return rows
            rows.append(row)

    def finish(self):
        """Require the source terminator and reject additional active rows after the supported configuration.

        Notes
        -----
        Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
        """
        self.expect([999])
        if self.i != len(self.rows):
            raise ValueError(f'{self.path}: unread input after end marker')


def ints(values):
    """Convert numeric values to integers only when every value is exactly integral.

    Parameters
    ----------
    values
        Numeric source values whose required integer/setting constraints are checked.

    Returns
    -------
    List of exact integers; fractional input is rejected rather than truncated.

    Notes
    -----
    Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
    """
    if any(v != int(v) for v in values):
        raise ValueError('Expected integer identifiers')
    return list(map(int, values))


def parameter(r, name, extended=True):
    """Read one control parameter record, validate its width, finite bounds and initialization, and retain its canonical name and source location.

    Parameters
    ----------
    r
        Stateful strict input Reader, advanced as parameter rows are consumed.
    name
        Exact section, output, field or stage identifier; matching is case sensitive.
    extended
        Whether the source parameter row uses the supported extended control layout.

    Returns
    -------
    Canonical parameter record with numeric columns and source provenance.

    Notes
    -----
    Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
    """
    row = r.take(14 if extended else 7)
    values = list(map(float, row.values))
    if not all(math.isfinite(v) for v in values):
        raise ValueError(f'{r.path}:{row.line}: nonfinite parameter')
    if not values[0] <= values[2] <= values[1]:
        raise ValueError(f'{r.path}:{row.line}: parameter initial outside bounds')
    return dict(name=name, columns=values, source_line=row.line, source_comment=row.comment)


def parse_data(path):
    """Parse the supported swordfish data layout, including original catches, indices, included/excluded length samples and calendar metadata.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.

    Returns
    -------
    Validated supported data dictionary with original observations.

    Notes
    -----
    Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
    """
    r = Reader(path)
    d = dict(start_year=r.integer(), end_year=r.integer(), seasons=r.integer())
    d['season_months'] = r.numbers(d['seasons'])
    d.update(subseasons=r.integer(), spawn_month=r.scalar(), sexes=r.integer(),
             plus_group_age=r.integer(), areas=r.integer(), nfleets=r.integer())
    if (d['sexes'], d['areas'], d['seasons'], d['nfleets'], d['plus_group_age']) != (2, 1, 4, 27, 15):
        raise ValueError('Reader is scoped to the supplied two-sex, one-area, quarterly model')
    if d['start_year'] > d['end_year'] or d['season_months'] != [3,3,3,3]:
        raise ValueError('Invalid time dimensions')
    fleets = []
    for f in range(1, d['nfleets'] + 1):
        row = r.take(6)
        spec = ints(list(map(float, row.values[:5])))
        fleets.append(dict(zip(['id','type','timing','area','catch_units','catch_multiplier','name'],
                               [f, *spec, row.values[5]])))
    catches = []
    for row in r.table(5):
        year, season, fleet = ints(row[:3])
        catches.append(dict(year=year, season=season, fleet=fleet, value=row[3], log_se=row[4],
                            initialization=(year == -999)))
    index_specs = [ints(r.numbers(4)) for _ in fleets]
    if [v[0] for v in index_specs] != list(range(1, d['nfleets'] + 1)):
        raise ValueError('Index definitions must follow numeric fleet order')
    indices = []
    for row in r.table(5):
        year, month, fleet, value, se = row
        ints([year,fleet])
        indices.append(dict(year=int(year), month=month, fleet=int(fleet), value=value, log_se=se,
                            use=(year > 0)))
    r.expect([0])  # no discard fleets
    r.expect([0])  # no mean body size observations
    r.expect([2])  # population bin grid specified by width and endpoints
    width, lo, hi = r.scalar(), r.scalar(), r.scalar()
    if width <= 0 or hi <= lo or (hi - lo) / width != int((hi - lo) / width):
        raise ValueError('Invalid population length grid')
    population_bins = [lo + i * width for i in range(int((hi - lo) / width) + 1)]
    r.expect([1])  # length compositions present
    comp_controls = [r.numbers(7) for _ in fleets]
    nbins = r.integer()
    bins = r.numbers(nbins)
    if bins != sorted(set(bins)):
        raise ValueError('Length bins must strictly increase')
    if bins != [10+5*i for i in range(51)] or population_bins != bins:
        raise ValueError('Unsupported length grid: expected common 51-bin 10-260 cm grid')
    compositions = []
    for row in r.table(6 + d['sexes'] * nbins):
        year, month, fleet, sex, part, sample_size = row[:6]
        ints([year, fleet, sex, part])
        compositions.append(dict(year=int(year), month=month, fleet=int(fleet), sex=int(sex),
                                 partition=int(part), sample_size=sample_size,
                                 use=(year > 0), values=row[6:]))
    # Explicitly reject additional observation types rather than ignoring them.
    for _ in range(7):
        r.expect([0])  # ages, size-at-age, env, sizefreq, tags, morphcomp, selpriors
    r.finish()
    result = dict(dimensions=d, fleets=fleets, catches=catches, index_specs=index_specs,
                  indices=indices, population_length_bins=population_bins, length_bins=bins,
                  length_controls=comp_controls, length_compositions=compositions)
    validate_data(result)
    return result


def validate_data(data):
    """Check observation years, fleets, positive errors/sample sizes, finite counts and duplicate sample keys before generating C++ data.

    Parameters
    ----------
    data
        Parsed observation/dimension/fleet dictionary used to validate compatible controls.

    Notes
    -----
    Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
    """
    d = data['dimensions']
    fs = {f['id']: f for f in data['fleets']}
    keys = set()
    for c in data['catches']:
        if not all(math.isfinite(c[k]) for k in ['value','log_se']):
            raise ValueError('Nonfinite catch')
        if c['fleet'] not in fs or fs[c['fleet']]['type'] != 1:
            raise ValueError('Catch references a non-catch fleet')
        if not 1 <= c['season'] <= d['seasons'] or c['value'] < 0 or c['log_se'] < 0:
            raise ValueError('Invalid catch observation')
        if not c['initialization'] and not d['start_year'] <= c['year'] <= d['end_year']:
            raise ValueError('Catch year outside model')
        key = (c['year'], c['season'], c['fleet'])
        if key in keys:
            raise ValueError('Duplicate catch observation')
        keys.add(key)
    for kind in ['indices', 'length_compositions']:
        keys = set()
        for x in data[kind]:
            if not d['start_year'] <= abs(x['year']) <= d['end_year']:
                raise ValueError(f'{kind}: year outside model')
            if not math.isfinite(x['month']):
                raise ValueError('Nonfinite observation month')
            if x['fleet'] not in fs or not 1 <= x['month'] < 13:
                raise ValueError(f'{kind}: invalid fleet/month')
            if kind == 'indices':
                if not all(math.isfinite(x[k]) for k in ['value','log_se']) or x['value'] <= 0 or x['log_se'] <= 0:
                    raise ValueError('Invalid lognormal index')
                key = (abs(x['year']), x['month'], x['fleet'])
            else:
                if x['sex'] not in (0,1,2,3) or x['partition'] not in (0,1,2):
                    raise ValueError('Invalid composition sex/partition code')
                if not math.isfinite(x['sample_size']) or x['sample_size'] <= 0 or any(not math.isfinite(v) or v < 0 for v in x['values']) or sum(x['values']) <= 0:
                    raise ValueError('Invalid composition counts/sample size')
                key = (abs(x['year']),x['month'],x['fleet'],x['sex'],x['partition'])
            if key in keys:
                raise ValueError(f'{kind}: duplicate observation')
            keys.add(key)


BIO_NAMES = ['L_at_Amin_Fem_GP_1','L_at_Amax_Fem_GP_1','VonBert_K_Fem_GP_1',
             'CV_young_Fem_GP_1','CV_old_Fem_GP_1','Wtlen_1_Fem_GP_1','Wtlen_2_Fem_GP_1',
             'Mat50%_Fem_GP_1','Mat_slope_Fem_GP_1','Eggs/kg_inter_Fem_GP_1',
             'Eggs/kg_slope_wt_Fem_GP_1','L_at_Amin_Mal_GP_1','L_at_Amax_Mal_GP_1',
             'VonBert_K_Mal_GP_1','CV_young_Mal_GP_1','CV_old_Mal_GP_1',
             'Wtlen_1_Mal_GP_1','Wtlen_2_Mal_GP_1','CohortGrowDev','FracFemale_GP_1']
SR_NAMES = ['SR_LN(R0)','SR_BH_steep','SR_sigmaR','SR_regime','SR_autocorr']


def parse_control(path, data):
    """Parse supported biological, selectivity, recruitment, catchability, hybrid-F and likelihood settings; reject unimplemented flags/layouts.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.
    data
        Parsed observation/dimension/fleet dictionary used to validate compatible controls.

    Returns
    -------
    Supported control dictionary compatible with the supplied data schema.

    Notes
    -----
    Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
    """
    r = Reader(path)
    d = data['dimensions']
    for expected in [0,1,1,4,1,1,0]:
        r.expect([expected])
    settlement = r.numbers(4)
    if settlement != [1,7,1,0]:
        raise ValueError('Unsupported settlement assignment')
    nblocks = r.integer()
    block_counts = ints(r.numbers(nblocks))
    blocks = [r.numbers(2*n) for n in block_counts]
    r.expect([1]); r.expect([0,0,0,0,1]); r.expect([3])
    mortality = [r.numbers(d['plus_group_age'] + 1) for _ in range(d['sexes'])]
    r.expect([1])  # von Bertalanffy growth
    a1, a2, plus_growth = r.scalar(), r.scalar(), r.scalar()
    r.expect([0]); r.expect([0]); r.expect([0]); r.expect([1])
    mature_age = r.integer()
    r.expect([1]); r.expect([0]); r.expect([1])
    parameters = [parameter(r, name) for name in BIO_NAMES]
    r.expect([0]*10)
    r.expect([3]); r.expect([1]); r.expect([1])
    parameters += [parameter(r, name) for name in SR_NAMES]
    r.expect([1])  # deviation vector recruitment model
    recruitment = dict(start=r.integer(), end=r.integer(), phase=r.integer())
    r.expect([1])
    recruitment['advanced'] = [r.scalar() for _ in range(13)]
    ballpark, ballpark_year = r.scalar(), r.integer()
    r.expect([3])
    max_f, hybrid_iterations = r.scalar(), r.integer()
    # One free initial-F parameter is associated with positive equilibrium catches.
    initial_pairs = [(x['fleet'],x['season']) for x in data['catches']
                     if x['initialization'] and x['value'] > 0]
    for f,s in initial_pairs:
        parameters.append(parameter(r,f'InitF_seas_{s}_flt_{f}',False))
    q_setup = r.table(6)
    for q in q_setup:
        parameters.append(parameter(r,f'LnQ_fleet_{int(q[0])}'))
    size_patterns = [ints(r.numbers(4)) for _ in data['fleets']]
    age_patterns = [ints(r.numbers(4)) for _ in data['fleets']]
    if any(row != [0,0,0,0] for row in age_patterns):
        raise ValueError('Nontrivial age selectivity is not supported by this reader')
    size_params = []
    for f, row in enumerate(size_patterns,1):
        pattern, discard, male, special = row
        if discard or male or pattern not in (0,1,15,24):
            raise ValueError(f'Unsupported size selectivity: fleet {f}, {row}')
        if pattern == 15 and not 1 <= special < f:
            raise ValueError('Mirror must reference an earlier fleet')
        n = {0:0,1:2,15:0,24:6}[pattern]
        group = [parameter(r,f'SizeSel_fleet_{f}_p{p}') for p in range(1,n+1)]
        parameters += group
        size_params += group
    # Expanded, manually supplied time-block parameters: read in structural order.
    for base in size_params:
        spec = base['columns']
        block = int(spec[12])
        if block:
            if not 1 <= block <= nblocks:
                raise ValueError('Invalid block assignment')
            for b in range(block_counts[block-1]):
                parameters.append(parameter(r,f"{base['name']}_block_{block}_{b+1}",False))
    r.expect([0]); r.expect([0])  # 2D AR1 and custom tagging disabled
    variance = r.table(3)
    maxphase, sd_offset = r.integer(), r.integer()
    lambdas = r.table(5)
    r.expect([0])  # no extra SD output definitions
    r.finish()
    return dict(settlement=settlement, blocks=blocks, natural_mortality=mortality,
                growth_reference_ages=[a1,a2], plus_group_growth=plus_growth,
                first_mature_age=mature_age, parameters=parameters, recruitment=recruitment,
                ballpark_f=ballpark, ballpark_year=ballpark_year, f_method=3,
                max_f=max_f, hybrid_iterations=hybrid_iterations, q_setup=q_setup,
                size_patterns=size_patterns, age_patterns=age_patterns,
                variance_adjustments=variance, max_lambda_phase=maxphase, sd_offset=sd_offset,
                lambda_changes=lambdas)


def load_inputs(raw):
    """Read and validate the paired data/control files, retaining hashes and explicit control/report distinctions for preparation.

    Parameters
    ----------
    raw
        Directory containing the eight original assessment inputs and reference report.

    Returns
    -------
    Combined original data/control/provenance dictionary.

    Notes
    -----
    Only the supplied SS3 v3.30.20 layout and explicitly implemented settings are accepted. Source comments are retained as provenance, never interpreted as executable instructions. Native control initial values and fitted Report.sso values remain distinct.
    """
    raw = Path(raw)
    data = parse_data(raw/'swo2023_v004.dat')
    control = parse_control(raw/'swo2023_v007.ctl',data)
    sources = {name: hashlib.sha256((raw/name).read_bytes()).hexdigest()
               for name in ['swo2023_v004.dat','swo2023_v007.ctl']}
    return dict(schema_version=1, source_sha256=sources, data=data, control=control,
                units=dict(numbers='thousands of fish', weight='kg', biomass='metric tonnes',
                           length='cm EFL', mortality='annual instantaneous rate'),
                status='Parsed original observations/control; no fit or likelihood evaluation')
