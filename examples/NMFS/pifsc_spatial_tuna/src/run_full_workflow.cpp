#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <process.h>
#else
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace {

#ifdef _WIN32
// The UCRT spawn functions rebuild a Windows command line from argv. Quote
// arguments explicitly so paths such as "2026 review/2026 quadra" survive
// the Windows C runtime's command-line parsing unchanged. Backslashes before
// a quote (and trailing backslashes before the closing quote) require doubling.
std::string quote_windows_argument(const std::string &argument) {
  if (argument.empty())
    return "\"\"";
  if (argument.find_first_of(" \t\n\v\"") == std::string::npos)
    return argument;

  std::string quoted = "\"";
  std::size_t backslashes = 0;
  for (char character : argument) {
    if (character == '\\') {
      ++backslashes;
      continue;
    }
    if (character == '"') {
      quoted.append(backslashes * 2 + 1, '\\');
      quoted.push_back('"');
    } else {
      quoted.append(backslashes, '\\');
      quoted.push_back(character);
    }
    backslashes = 0;
  }
  quoted.append(backslashes * 2, '\\');
  quoted.push_back('"');
  return quoted;
}
#endif

struct WorkflowOptions {
  std::string config = "config/tuna_assessment.conf";
  std::string output_dir = "build/assessment_outputs/data";
  std::string report_dir = "build/assessment_outputs/report";
  std::string assessment_binary;
  std::string report_script = "scripts/build_assessment_report.R";
  std::string rscript = "Rscript";
  int simulations = 50;
};

std::string require_value(int &index, int argc, char **argv,
                          const std::string &option) {
  if (++index >= argc)
    throw std::invalid_argument(option + " requires a value");
  return argv[index];
}

WorkflowOptions parse_options(int argc, char **argv) {
  WorkflowOptions options;
  const std::filesystem::path self =
      std::filesystem::absolute(argv[0]).lexically_normal();
  options.assessment_binary =
      (self.parent_path() / "advanced_tuna_spatial_assessment_example")
          .string();
#ifdef _WIN32
  options.assessment_binary += ".exe";
#endif

  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--config")
      options.config = require_value(i, argc, argv, argument);
    else if (argument == "--output-dir")
      options.output_dir = require_value(i, argc, argv, argument);
    else if (argument == "--report-dir")
      options.report_dir = require_value(i, argc, argv, argument);
    else if (argument == "--assessment-binary")
      options.assessment_binary = require_value(i, argc, argv, argument);
    else if (argument == "--report-script")
      options.report_script = require_value(i, argc, argv, argument);
    else if (argument == "--rscript")
      options.rscript = require_value(i, argc, argv, argument);
    else if (argument == "--simulations") {
      options.simulations = std::stoi(require_value(i, argc, argv, argument));
      if (options.simulations <= 0)
        throw std::invalid_argument("--simulations must be positive");
    } else if (argument == "--help" || argument == "-h") {
      std::cout
          << "Usage: run_tuna_full_workflow [options]\n\n"
          << "  --config PATH              Assessment configuration\n"
          << "  --output-dir PATH          Machine-readable output directory\n"
          << "  --report-dir PATH          Generated report directory\n"
          << "  --simulations N            Simulation-estimation replicates\n"
          << "  --assessment-binary PATH   Assessment executable override\n"
          << "  --report-script PATH       Report-builder script override\n"
          << "  --rscript PATH             Rscript executable override\n";
      std::exit(0);
    } else {
      throw std::invalid_argument("unknown option: " + argument);
    }
  }
  return options;
}

int run_process(
    const std::string &stage, const std::vector<std::string> &args,
    const std::vector<std::pair<std::string, std::string>> &environment = {}) {
  std::cout << "\n[workflow]\n  stage: " << stage << "\n  status: running\n"
            << std::flush;
#ifdef _WIN32
  struct SavedEnvironment {
    std::string name;
    std::string value;
    bool existed = false;
  };
  std::vector<SavedEnvironment> saved_environment;
  saved_environment.reserve(environment.size());
  for (const auto &entry : environment) {
    const char *old_value = std::getenv(entry.first.c_str());
    saved_environment.push_back({entry.first,
                                 old_value == nullptr ? "" : old_value,
                                 old_value != nullptr});
    if (_putenv_s(entry.first.c_str(), entry.second.c_str()) != 0)
      throw std::runtime_error("could not set environment variable " +
                               entry.first + " for " + stage);
  }

  std::vector<std::string> quoted_args;
  quoted_args.reserve(args.size());
  for (const std::string &argument : args)
    quoted_args.push_back(quote_windows_argument(argument));

  std::vector<const char *> child_argv;
  child_argv.reserve(args.size() + 1);
  for (const std::string &argument : quoted_args)
    child_argv.push_back(argument.c_str());
  child_argv.push_back(nullptr);

  const intptr_t process_status =
      _spawnvp(_P_WAIT, args[0].c_str(), child_argv.data());
  const int spawn_errno = errno;

  for (const auto &entry : saved_environment) {
    _putenv_s(entry.name.c_str(), entry.existed ? entry.value.c_str() : "");
  }

  if (process_status == -1)
    throw std::runtime_error("could not launch " + stage +
                             ": errno=" + std::to_string(spawn_errno));
  const int exit_code = static_cast<int>(process_status);
#else
  const pid_t child = fork();
  if (child < 0)
    throw std::runtime_error("fork failed for " + stage +
                             ": errno=" + std::to_string(errno));
  if (child == 0) {
    for (const auto &entry : environment) {
      if (setenv(entry.first.c_str(), entry.second.c_str(), 1) != 0)
        _exit(126);
    }
    std::vector<char *> child_argv;
    child_argv.reserve(args.size() + 1);
    for (const std::string &argument : args)
      child_argv.push_back(const_cast<char *>(argument.c_str()));
    child_argv.push_back(nullptr);
    execvp(child_argv[0], child_argv.data());
    _exit(errno == ENOENT ? 127 : 126);
  }

  int status = 0;
  while (waitpid(child, &status, 0) < 0) {
    if (errno != EINTR)
      throw std::runtime_error("waitpid failed for " + stage +
                               ": errno=" + std::to_string(errno));
  }
  const int exit_code =
      WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
#endif
  std::cout << "\n[workflow]\n  stage: " << stage
            << "\n  status: " << (exit_code == 0 ? "complete" : "failed")
            << "\n  exit_code: " << exit_code << "\n"
            << std::flush;
  return exit_code;
}

} // namespace

int main(int argc, char **argv) {
  try {
    const WorkflowOptions options = parse_options(argc, argv);
    std::filesystem::create_directories(options.output_dir);
    std::filesystem::create_directories(options.report_dir);

    int status = run_process(
        "fit_and_simulation",
        {options.assessment_binary, "--config", options.config},
        {{"QUADRA_TUNA_BASELINE_ONLY", "0"},
         {"QUADRA_TUNA_OUTPUT_DIR", options.output_dir},
         {"QUADRA_TUNA_SIMULATIONS", std::to_string(options.simulations)}});
    if (status != 0)
      return status;

    status = run_process("report_generation",
                         {options.rscript, options.report_script,
                          options.output_dir, options.report_dir});
    if (status != 0)
      return status;

    std::cout << "\n[workflow]\n"
              << "  status: complete\n"
              << "  data_dir: " << options.output_dir << "\n"
              << "  report_dir: " << options.report_dir << "\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "full tuna workflow failed: " << error.what() << '\n';
    return 2;
  }
}
