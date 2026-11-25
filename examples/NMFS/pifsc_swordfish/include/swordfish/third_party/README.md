# Header-only runtime dependencies

These pinned sources are bundled so building/running needs no package install or
network access. Neither dependency launches external processes.

- `json.hpp`: JSON for Modern C++ **3.11.3**, from
  https://github.com/nlohmann/json/blob/v3.11.3/single_include/nlohmann/json.hpp
  (MIT; see `JSON_LICENSE.MIT`). Used for strict native JSON intake and output.
- `picosha2.h`: PicoSHA2 **1.0.0**, from
  https://github.com/okdshin/PicoSHA2/blob/v1.0.0/picosha2.h
  (MIT; full license retained in the header). Used for provenance SHA-256 hashes.
