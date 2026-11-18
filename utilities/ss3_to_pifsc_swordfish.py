#!/usr/bin/env python3
"""Convert the supported NP swordfish SS3 inputs into native Quadra JSON."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       'examples/NMFS/pifsc_swordfish/scripts'))
from ss3_to_native import main

if __name__ == '__main__':
    main()
