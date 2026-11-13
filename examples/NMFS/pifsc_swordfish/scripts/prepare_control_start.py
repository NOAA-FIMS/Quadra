"""Generate the cold-start factory using only original .dat/.ctl inputs, never Report.sso."""
from pathlib import Path
from swordfish_inputs import load_inputs
from prepare_model import generate_header


def main():
    """Validate the two source inputs and emit a dedicated original-control C++ factory.

    The ordinary audit preparation also extracts Report.sso fixtures. This entry
    point deliberately omits that path so cold-start builds need no SS3 fitted
    output or covariance file at either preparation or execution time.
    """
    root=Path(__file__).resolve().parents[1]
    output=root/'build';output.mkdir(exist_ok=True)
    generate_header(load_inputs(root/'data/raw'),output/'swordfish_control_input.hpp',factory_name='imported_swordfish_control_data')
    print('Control-only input factory generated from .dat/.ctl (no Report.sso).')


if __name__=='__main__':main()
