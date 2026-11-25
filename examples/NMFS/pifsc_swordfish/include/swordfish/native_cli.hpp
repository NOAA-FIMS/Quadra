/** @file native_cli.hpp Standalone conversion, evaluation and fitting options.
 */
#pragma once
#include "csv.hpp"
#include "native_files.hpp"
#include <iostream>
#include <map>
#include <set>
namespace swordfish {
struct CLIError : std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};
struct NativeOptions {
  bool help = false, convert = false, legacy = false, uncertainty = false;
  fs::path input, output, data, control, parameters, metric;
  std::string mode = "fit";
  int iterations = 200;
  double tolerance = 1e-4;
};
inline NativeOptions native_options(int argc, char **argv) {
  NativeOptions o;
  if (argc == 1) {
    o.help = true;
    return o;
  }
  if (argv[1][0] != '-') {
    o.legacy = true;
    if (argc < 5 || argc > 7)
      throw CLIError(
          "Backend usage: OUTPUT MODE MAX_ITERATIONS TOLERANCE [uncertainty] "
          "[metric=PATH]. Run --help for the full menu.");
    o.output = argv[1];
    o.mode = argv[2];
    try {
      o.iterations = csv_integer(argv[3]);
      o.tolerance = csv_number(argv[4]);
    } catch (const std::exception &) {
      throw CLIError("Invalid optimizer limits");
    }
    for (int i = 5; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "uncertainty" && !o.uncertainty)
        o.uncertainty = true;
      else if (arg.rfind("metric=", 0) == 0 && o.metric.empty() &&
               arg.size() > 7)
        o.metric = arg.substr(7);
      else
        throw CLIError("Unknown/duplicate backend option: " + arg);
    }
  } else {
    o.convert = std::string(argv[1]) == "--convert-ss3";
    std::map<std::string, std::string> values;
    std::set<std::string> seen;
    const std::set<std::string> options =
        o.convert ? std::set<std::string>{"--data", "--control", "--parameters",
                                          "--output"}
                  : std::set<std::string>{"--input",
                                          "--output",
                                          "--mode",
                                          "--max-iterations",
                                          "--gradient-tolerance",
                                          "--optimizer-metric"};
    for (int i = o.convert ? 2 : 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--help" || arg == "-h") {
        o.help = true;
        continue;
      }
      if (!seen.insert(arg).second)
        throw CLIError("Duplicate option: " + arg);
      if (arg == "--uncertainty" && !o.convert) {
        o.uncertainty = true;
        continue;
      }
      if (!options.count(arg))
        throw CLIError("Unknown option: " + arg + ". Run --help.");
      if (i + 1 == argc || std::string(argv[i + 1]).rfind("--", 0) == 0)
        throw CLIError("Missing value for " + arg);
      values[arg] = argv[++i];
      if (values[arg].empty())
        throw CLIError("Empty value for " + arg);
    }
    if (o.help)
      return o;
    auto required = [&](const std::string &name) {
      if (!values.count(name))
        throw CLIError("Required option: " + name + ". Run --help.");
      return values.at(name);
    };
    if (o.convert) {
      o.data = required("--data");
      o.control = required("--control");
      o.parameters = values["--parameters"];
    } else {
      o.input = required("--input");
      if (values.count("--mode"))
        o.mode = values.at("--mode");
      try {
        if (values.count("--max-iterations"))
          o.iterations = csv_integer(values.at("--max-iterations"));
        if (values.count("--gradient-tolerance"))
          o.tolerance = csv_number(values.at("--gradient-tolerance"));
      } catch (const std::exception &) {
        throw CLIError("Invalid optimizer limits");
      }
      o.metric = values["--optimizer-metric"];
    }
    o.output = required("--output");
  }
  if (!o.convert) {
    if (o.mode != "fit" && o.mode != "evaluate")
      throw CLIError("--mode must be fit or evaluate");
    if (o.mode != "fit" && (o.uncertainty || !o.metric.empty()))
      throw CLIError("--uncertainty and --optimizer-metric require --mode fit");
    if (o.iterations < 1 || !std::isfinite(o.tolerance) || !(o.tolerance > 0))
      throw CLIError(
          "Optimizer iterations and tolerance must be positive and finite");
  }
  return o;
}
inline void print_native_help(std::ostream &out, bool conversion = false) {
  out << R"HELP(pifsc_swordfish - Quadra

USAGE
  pifsc_swordfish --input MODEL.json --output NEW_DIRECTORY [model options]
  pifsc_swordfish --convert-ss3 --data FILE --control FILE [--parameters FILE]
                   --output NEW.json
  pifsc_swordfish --help | -h

MODEL OPTIONS
  --input FILE              Native JSON model (quadra.swordfish.v1; required).
  --output DIRECTORY        New output directory (required; never overwritten).
  --mode fit|evaluate       Fit the model (default) or evaluate its starts.
  --max-iterations N        Optimizer limit (default: 200; positive integer).
  --gradient-tolerance X    Final gradient tolerance (default: 0.0001).
  --uncertainty             Local reference-point intervals (requires fit).
  --optimizer-metric FILE   Labeled covariance for optimizer scaling (fit).

CONVERT SS3 FILES TO NATIVE INPUT
  --convert-ss3             Select conversion mode; place this option first.
  --data FILE               SS3 observation (.dat) file (required).
  --control FILE            SS3 control (.ctl) file (required).
  --parameters FILE         Optional SS3 .par file supplying initial values.
  --output FILE             New native JSON file (required).

  Writes JSON plus a .provenance.json sidecar; existing files are protected.
  Supports the North Pacific swordfish SS3 3.30.20 layout.
  Without --parameters, literal control starts and zero deviations are kept;
  those starts require the separate catch-feasibility workflow before fitting.

EXAMPLES (from the pifsc_swordfish directory)
  ./build/pifsc_swordfish --convert-ss3 --data data/raw/swo2023_v004.dat \
    --control data/raw/swo2023_v007.ctl --parameters data/raw/ss.par \
    --output build/native/from-ss3.json

  ./build/pifsc_swordfish --input build/native/from-ss3.json \
    --output build/native/evaluation --mode evaluate

  ./build/pifsc_swordfish --input build/native/from-ss3.json \
    --output build/native/my-fit --uncertainty

OUTPUT AND PROVENANCE
  Runs save model_input.json (exact input bytes), run_status.json with SHA-256
  input/binary/metric hashes, fitted outputs and derived reference points.
  Conversion, JSON validation, fitting and reporting run entirely in C++.
  No Python installation, converter scripts or environment settings are needed.

LEGACY BACKEND MODE (compatibility with existing scripts)
  pifsc_swordfish OUTPUT MODE MAX_ITERATIONS TOLERANCE [uncertainty] [metric=PATH]
  Reads the private typed stream from stdin. Named --input reads JSON directly.

MORE HELP
  pifsc_swordfish --convert-ss3 --help
  See docs/NATIVE_INPUT.md for schema, initialization, and output details.
)HELP";
  (void)conversion;
}
} // namespace swordfish
