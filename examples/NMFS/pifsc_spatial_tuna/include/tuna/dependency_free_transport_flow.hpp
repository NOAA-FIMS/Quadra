#pragma once

#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace quadra::transport {
struct DenseLayer {
  std::size_t input_m = 0;
  std::size_t output_m = 0;
  std::vector<float> weight_m;
  std::vector<float> bias_m;

  std::vector<double> operator()(const std::vector<double> &input) const {
    std::vector<double> output(output_m);
    for (std::size_t row = 0; row < output_m; ++row) {
      double value = bias_m[row];
      for (std::size_t column = 0; column < input_m; ++column)
        value += weight_m[row * input_m + column] * input[column];
      output[row] = value;
    }
    return output;
  }
};

struct CouplingWeights {
  DenseLayer input_m;
  DenseLayer middle_m;
  DenseLayer output_m;

  std::vector<double> operator()(const std::vector<double> &x) const {
    std::vector<double> hidden = input_m(x);
    for (double &value : hidden)
      value = std::max(0.0, value);
    hidden = middle_m(hidden);
    for (double &value : hidden)
      value = std::max(0.0, value);
    return output_m(hidden);
  }
};

class DependencyFreeRealNVP {
public:
  void load(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    char magic[8]{};
    input.read(magic, sizeof(magic));
    if (!input || std::string(magic, sizeof(magic)) != "QFLOW001")
      throw std::runtime_error("invalid dependency-free flow archive: " + path);
    const std::uint32_t version = read<std::uint32_t>(input);
    if (version != 1)
      throw std::runtime_error("unsupported QFLOW version");
    dimension_m = read<std::uint32_t>(input);
    const std::uint32_t layers = read<std::uint32_t>(input);
    hidden_m = read<std::uint32_t>(input);
    log_scale_limit_m = read<double>(input);
    if (dimension_m < 2 || layers < 2 || hidden_m < 2)
      throw std::runtime_error("invalid QFLOW architecture");
    mean_m = read_vector(input, dimension_m);
    scale_m = read_vector(input, dimension_m);
    networks_m.clear();
    networks_m.reserve(layers);
    for (std::uint32_t layer = 0; layer < layers; ++layer) {
      CouplingWeights network;
      network.input_m = read_layer(input, dimension_m, hidden_m);
      network.middle_m = read_layer(input, hidden_m, hidden_m);
      network.output_m = read_layer(input, hidden_m, 2 * dimension_m);
      networks_m.push_back(std::move(network));
    }
    if (input.peek() != std::ifstream::traits_type::eof())
      throw std::runtime_error("unexpected trailing data in QFLOW archive");
  }

  std::size_t dimension() const { return dimension_m; }
  std::size_t layers() const { return networks_m.size(); }
  std::size_t hidden() const { return hidden_m; }

  std::vector<double> forward_noise(const std::vector<double> &input) const {
    if (input.size() != dimension_m)
      throw std::invalid_argument("QFLOW forward dimension mismatch");
    std::vector<double> x(dimension_m);
    for (std::size_t j = 0; j < dimension_m; ++j)
      x[j] = (input[j] - mean_m[j]) / scale_m[j];
    for (std::size_t layer = 0; layer < networks_m.size(); ++layer) {
      std::vector<double> masked(dimension_m);
      for (std::size_t j = 0; j < dimension_m; ++j)
        if ((j + layer) % 2 == 0)
          masked[j] = x[j];
      const auto parameters = networks_m[layer](masked);
      for (std::size_t j = 0; j < dimension_m; ++j)
        if ((j + layer) % 2 != 0) {
          const double s = log_scale_limit_m * std::tanh(parameters[j]);
          x[j] = x[j] * std::exp(s) + parameters[dimension_m + j];
        }
    }
    return x;
  }

  double log_density(const std::vector<double> &input) const {
    if (input.size() != dimension_m)
      throw std::invalid_argument("QFLOW density dimension mismatch");
    std::vector<double> x(dimension_m);
    double log_jacobian = 0.0;
    for (std::size_t j = 0; j < dimension_m; ++j) {
      x[j] = (input[j] - mean_m[j]) / scale_m[j];
      log_jacobian -= std::log(scale_m[j]);
    }
    for (std::size_t layer = 0; layer < networks_m.size(); ++layer) {
      std::vector<double> masked(dimension_m);
      for (std::size_t j = 0; j < dimension_m; ++j)
        if ((j + layer) % 2 == 0)
          masked[j] = x[j];
      const std::vector<double> parameters = networks_m[layer](masked);
      for (std::size_t j = 0; j < dimension_m; ++j)
        if ((j + layer) % 2 != 0) {
          const double log_scale = log_scale_limit_m * std::tanh(parameters[j]);
          x[j] = x[j] * std::exp(log_scale) + parameters[dimension_m + j];
          log_jacobian += log_scale;
        }
    }
    constexpr double log_two_pi = 1.83787706640934548356;
    double base = 0.0;
    for (const double value : x)
      base -= 0.5 * (value * value + log_two_pi);
    return base + log_jacobian;
  }

  std::vector<double> inverse(const std::vector<double> &noise) const {
    if (noise.size() != dimension_m)
      throw std::invalid_argument("QFLOW inverse dimension mismatch");
    std::vector<double> x = noise;
    for (std::size_t reverse = networks_m.size(); reverse-- > 0;) {
      std::vector<double> masked(dimension_m);
      for (std::size_t j = 0; j < dimension_m; ++j)
        if ((j + reverse) % 2 == 0)
          masked[j] = x[j];
      const std::vector<double> parameters = networks_m[reverse](masked);
      for (std::size_t j = 0; j < dimension_m; ++j)
        if ((j + reverse) % 2 != 0) {
          const double log_scale = log_scale_limit_m * std::tanh(parameters[j]);
          x[j] = (x[j] - parameters[dimension_m + j]) * std::exp(-log_scale);
        }
    }
    for (std::size_t j = 0; j < dimension_m; ++j)
      x[j] = mean_m[j] + scale_m[j] * x[j];
    return x;
  }

private:
  template <class Type> static Type read(std::ifstream &input) {
    Type value{};
    input.read(reinterpret_cast<char *>(&value), sizeof(value));
    if (!input)
      throw std::runtime_error("truncated QFLOW archive");
    return value;
  }

  static std::vector<float> read_vector(std::ifstream &input,
                                        std::size_t count) {
    std::vector<float> values(count);
    input.read(reinterpret_cast<char *>(values.data()),
               static_cast<std::streamsize>(count * sizeof(float)));
    if (!input)
      throw std::runtime_error("truncated QFLOW vector");
    for (const float value : values)
      if (!std::isfinite(value))
        throw std::runtime_error("non-finite QFLOW weight");
    return values;
  }

  static DenseLayer read_layer(std::ifstream &input, std::size_t columns,
                               std::size_t rows) {
    DenseLayer layer;
    layer.input_m = columns;
    layer.output_m = rows;
    layer.weight_m = read_vector(input, rows * columns);
    layer.bias_m = read_vector(input, rows);
    return layer;
  }

  std::size_t dimension_m = 0;
  std::size_t hidden_m = 0;
  double log_scale_limit_m = 1.5;
  std::vector<float> mean_m;
  std::vector<float> scale_m;
  std::vector<CouplingWeights> networks_m;
};
} // namespace quadra::transport
