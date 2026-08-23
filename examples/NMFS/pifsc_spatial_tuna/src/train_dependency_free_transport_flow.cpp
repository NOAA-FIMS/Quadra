#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "../include/tuna/dependency_free_transport_flow.hpp"

namespace {
using Matrix =
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Row = Eigen::Matrix<float, 1, Eigen::Dynamic, Eigen::RowMajor>;

struct Options {
  std::vector<std::string> draws;
  std::string output, geometry, fit;
  int layers = 8, hidden = 64, epochs = 1000, patience = 120;
  uint64_t seed = 20260823;
};
std::vector<std::string> split(const std::string &s, char c) {
  std::vector<std::string> o;
  std::istringstream in(s);
  std::string x;
  while (std::getline(in, x, c))
    o.push_back(x);
  return o;
}
Options options(int n, char **v) {
  Options o;
  for (int i = 1; i < n; ++i) {
    std::string a = v[i];
    auto next = [&]() {
      if (++i >= n)
        throw std::runtime_error("missing " + a);
      return std::string(v[i]);
    };
    if (a == "--draws")
      while (i + 1 < n && std::string(v[i + 1]).rfind("--", 0) != 0)
        o.draws.push_back(v[++i]);
    else if (a == "--output")
      o.output = next();
    else if (a == "--geometry-cache")
      o.geometry = next();
    else if (a == "--fit-metadata")
      o.fit = next();
    else if (a == "--layers")
      o.layers = std::stoi(next());
    else if (a == "--hidden")
      o.hidden = std::stoi(next());
    else if (a == "--epochs")
      o.epochs = std::stoi(next());
    else if (a == "--patience")
      o.patience = std::stoi(next());
    else if (a == "--seed")
      o.seed = std::stoull(next());
    else
      throw std::runtime_error("unknown argument " + a);
  }
  if (o.draws.empty() || o.output.empty() || o.geometry.empty() ||
      o.fit.empty())
    throw std::runtime_error("missing required arguments");
  return o;
}
std::string hash_file(const std::string &p) {
  std::ifstream in(p, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot hash " + p);
  uint64_t h = 1469598103934665603ULL;
  char b[65536];
  while (in) {
    in.read(b, sizeof b);
    for (std::streamsize i = 0; i < in.gcount(); ++i) {
      h ^= (unsigned char)b[i];
      h *= 1099511628211ULL;
    }
  }
  std::ostringstream s;
  s << std::hex << std::setw(16) << std::setfill('0') << h;
  return s.str();
}
std::string fingerprint(const std::string &p) {
  std::ifstream in(p);
  std::string a, b;
  std::getline(in, a);
  if (a.rfind("fingerprint,converged", 0) == 0) {
    std::getline(in, b);
    return split(b, ',').at(0);
  }
  do {
    auto f = split(a, ',');
    if (f.size() > 1 && f[0] == "fingerprint")
      return f[1];
  } while (std::getline(in, a));
  throw std::runtime_error("missing fingerprint");
}
struct Corpus {
  std::vector<std::string> names;
  Matrix x;
  std::vector<std::pair<int, int>> groups;
};
Corpus corpus(const std::vector<std::string> &paths) {
  using Key = std::tuple<int, int, int>;
  std::map<Key, std::unordered_map<std::string, float>> g;
  std::vector<std::string> names;
  std::set<std::string> known;
  for (size_t s = 0; s < paths.size(); ++s) {
    std::ifstream in(paths[s]);
    std::string line;
    std::getline(in, line);
    auto h = split(line, ',');
    auto pos = [&](std::string n) {
      return size_t(std::find(h.begin(), h.end(), n) - h.begin());
    };
    size_t cc = pos("chain"), ic = pos("iteration"), pc = pos("parameter"),
           vc = pos("value");
    while (std::getline(in, line)) {
      auto f = split(line, ',');
      g[{int(s), std::stoi(f[cc]), std::stoi(f[ic])}][f[pc]] = std::stof(f[vc]);
      if (known.insert(f[pc]).second)
        names.push_back(f[pc]);
    }
  }
  std::set<std::vector<float>> unique;
  std::vector<std::vector<float>> rows;
  std::vector<std::pair<int, int>> groups;
  for (auto &e : g) {
    std::vector<float> r;
    for (auto &n : names) {
      auto q = e.second.find(n);
      if (q == e.second.end())
        throw std::runtime_error("incomplete draw");
      r.push_back(q->second);
    }
    if (unique.insert(r).second) {
      rows.push_back(r);
      groups.push_back({std::get<0>(e.first), std::get<1>(e.first)});
    }
  }
  Matrix x(rows.size(), names.size());
  for (int i = 0; i < x.rows(); ++i)
    for (int j = 0; j < x.cols(); ++j)
      x(i, j) = rows[i][j];
  return {names, x, groups};
}

struct Parameter {
  Matrix x, m, v;
  Parameter() = default;
  Parameter(int r, int c)
      : x(r, c), m(Matrix::Zero(r, c)), v(Matrix::Zero(r, c)) {}
  void step(const Matrix &g, int t, float lr) {
    constexpr float b1 = .9f, b2 = .999f, eps = 1e-8f, wd = 1e-4f;
    m = b1 * m + (1 - b1) * g;
    v = b2 * v + (1 - b2) * g.array().square().matrix();
    float a = lr * std::sqrt(1 - std::pow(b2, t)) / (1 - std::pow(b1, t));
    x.array() *= (1 - lr * wd);
    x.array() -= a * m.array() / (v.array().sqrt() + eps);
  }
};
struct Net {
  Parameter w1, b1, w2, b2, w3, b3;
  Net() = default;
  Net(int d, int h)
      : w1(h, d), b1(1, h), w2(h, h), b2(1, h), w3(2 * d, h), b3(1, 2 * d) {}
};
struct Cache {
  std::vector<Matrix> x, z1, h1, z2, h2, p, s, e;
};
struct Gradients {
  Matrix w1, b1, w2, b2, w3, b3;
};
struct Flow {
  int d, h, l;
  float limit = 1.5f;
  Row mean, scale;
  std::vector<Net> n;
  Flow(int D, int L, int H, uint64_t seed)
      : d(D), h(H), l(L), mean(D), scale(D) {
    std::mt19937_64 rng(seed);
    n.reserve(l);
    for (int k = 0; k < l; ++k) {
      n.emplace_back(d, h);
      auto init = [&](Parameter &p, float sd) {
        std::normal_distribution<float> z(0, sd);
        for (int i = 0; i < p.x.size(); ++i)
          p.x.data()[i] = z(rng);
      };
      init(n.back().w1, std::sqrt(2.f / d));
      init(n.back().w2, std::sqrt(2.f / h));
      n.back().w3.x.setZero();
      n.back().b1.x.setZero();
      n.back().b2.x.setZero();
      n.back().b3.x.setZero();
    }
  }
  Row mask(int k) const {
    Row m(d);
    for (int j = 0; j < d; ++j)
      m[j] = ((j + k) % 2 == 0);
    return m;
  }
  float forward(const Matrix &raw, Cache *c = nullptr) const {
    Cache local;
    if (!c)
      c = &local;
    c->x.clear();
    c->z1.clear();
    c->h1.clear();
    c->z2.clear();
    c->h2.clear();
    c->p.clear();
    c->s.clear();
    c->e.clear();
    Matrix x = (raw.rowwise() - mean).array().rowwise() / scale.array();
    c->x.push_back(x);
    Matrix logj = Matrix::Constant(raw.rows(), 1, -scale.array().log().sum());
    for (int k = 0; k < l; ++k) {
      Row m = mask(k), u = Row::Ones(d) - m;
      Matrix masked = x.array().rowwise() * m.array();
      Matrix z1 = masked * n[k].w1.x.transpose();
      z1.rowwise() += n[k].b1.x.row(0);
      Matrix h1 = z1.cwiseMax(0);
      Matrix z2 = h1 * n[k].w2.x.transpose();
      z2.rowwise() += n[k].b2.x.row(0);
      Matrix h2 = z2.cwiseMax(0);
      Matrix p = h2 * n[k].w3.x.transpose();
      p.rowwise() += n[k].b3.x.row(0);
      Matrix s = (p.leftCols(d).array().tanh() * limit).matrix();
      s.array().rowwise() *= u.array();
      Matrix e = s.array().exp();
      Matrix t = p.rightCols(d).array().rowwise() * u.array();
      x = (x.array().rowwise() * m.array() +
           u.replicate(x.rows(), 1).array() *
               (x.array() * e.array() + t.array()))
              .matrix();
      logj += s.rowwise().sum();
      c->z1.push_back(z1);
      c->h1.push_back(h1);
      c->z2.push_back(z2);
      c->h2.push_back(h2);
      c->p.push_back(p);
      c->s.push_back(s);
      c->e.push_back(e);
      c->x.push_back(x);
    }
    return (.5f * x.array().square().rowwise().sum().mean() +
            .5f * d * 1.8378770664f - logj.mean());
  }
  std::vector<Gradients> gradients(const Matrix &raw) const {
    Cache c;
    forward(raw, &c);
    float inv = 1.f / raw.rows();
    Matrix G = c.x.back() * inv;
    std::vector<Gradients> gs(l);
    for (int k = l - 1; k >= 0; --k) {
      Row m = mask(k), u = Row::Ones(d) - m;
      Matrix gp = Matrix::Zero(raw.rows(), 2 * d);
      Matrix gscale = (G.array() * c.x[k].array() * c.e[k].array()).matrix();
      gscale.array() -= inv;
      gscale.array().rowwise() *= u.array();
      gp.leftCols(d) = (gscale.array() * limit *
                        (1 - c.p[k].leftCols(d).array().tanh().square()))
                           .matrix();
      gp.rightCols(d) = G.array().rowwise() * u.array();
      gs[k].w3 = gp.transpose() * c.h2[k];
      gs[k].b3 = gp.colwise().sum();
      Matrix gh2 = gp * n[k].w3.x;
      gh2.array() *= (c.z2[k].array() > 0).cast<float>();
      gs[k].w2 = gh2.transpose() * c.h1[k];
      gs[k].b2 = gh2.colwise().sum();
      Matrix gh1 = gh2 * n[k].w2.x;
      gh1.array() *= (c.z1[k].array() > 0).cast<float>();
      gs[k].w1 =
          gh1.transpose() * (c.x[k].array().rowwise() * m.array()).matrix();
      gs[k].b1 = gh1.colwise().sum();
      Matrix gin = gh1 * n[k].w1.x;
      G = (G.array().rowwise() * m.array() + gin.array().rowwise() * m.array() +
           G.array() * c.e[k].array() * u.replicate(G.rows(), 1).array())
              .matrix();
    }
    return gs;
  }
  void train_batch(const Matrix &raw, int step, float lr) {
    std::vector<Gradients> gs = gradients(raw);
    double norm = 0;
    for (auto &g : gs)
      for (Matrix *p : {&g.w1, &g.b1, &g.w2, &g.b2, &g.w3, &g.b3})
        norm += p->squaredNorm();
    float clip = norm > 100 ? float(10 / std::sqrt(norm)) : 1;
    for (int k = 0; k < l; ++k) {
      n[k].w1.step(gs[k].w1 * clip, step, lr);
      n[k].b1.step(gs[k].b1 * clip, step, lr);
      n[k].w2.step(gs[k].w2 * clip, step, lr);
      n[k].b2.step(gs[k].b2 * clip, step, lr);
      n[k].w3.step(gs[k].w3 * clip, step, lr);
      n[k].b3.step(gs[k].b3 * clip, step, lr);
    }
  }
};

double check_gradients(Flow &flow, const Matrix &batch) {
  const auto analytic = flow.gradients(batch);
  struct Probe {
    Parameter *parameter;
    const Matrix *gradient;
    int row;
    int column;
  };
  std::vector<Probe> probes = {
      {&flow.n.front().w3, &analytic.front().w3, 0, 0},
      {&flow.n.front().b3, &analytic.front().b3, 0, flow.d},
      {&flow.n.front().w2, &analytic.front().w2, 0, 0},
      {&flow.n.back().w1, &analytic.back().w1, 0, 0},
      {&flow.n.back().w3, &analytic.back().w3, flow.d, 0},
      {&flow.n.back().b3, &analytic.back().b3, 0, 0}};
  constexpr float step = 1e-3f;
  double maximum_scaled_error = 0.0;
  for (const auto &probe : probes) {
    float &value = probe.parameter->x(probe.row, probe.column);
    const float original = value;
    value = original + step;
    const double plus = flow.forward(batch);
    value = original - step;
    const double minus = flow.forward(batch);
    value = original;
    const double numerical = (plus - minus) / (2.0 * step);
    const double expected = (*probe.gradient)(probe.row, probe.column);
    maximum_scaled_error =
        std::max(maximum_scaled_error,
                 std::abs(numerical - expected) /
                     std::max(1.0, std::abs(numerical) + std::abs(expected)));
  }
  if (!std::isfinite(maximum_scaled_error) || maximum_scaled_error > 2e-2)
    throw std::runtime_error("RealNVP finite-difference gradient check failed");
  return maximum_scaled_error;
}

template <class T> void bin(std::ofstream &o, const T &v) {
  o.write((const char *)&v, sizeof v);
}
void matrix(std::ofstream &o, const Matrix &m) {
  o.write((const char *)m.data(), m.size() * sizeof(float));
}
void save(const Flow &f, const std::string &p) {
  std::ofstream o(p, std::ios::binary);
  o.write("QFLOW001", 8);
  bin(o, uint32_t(1));
  bin(o, uint32_t(f.d));
  bin(o, uint32_t(f.l));
  bin(o, uint32_t(f.h));
  bin(o, double(f.limit));
  o.write((const char *)f.mean.data(), f.d * 4);
  o.write((const char *)f.scale.data(), f.d * 4);
  for (auto &n : f.n) {
    matrix(o, n.w1.x);
    matrix(o, n.b1.x);
    matrix(o, n.w2.x);
    matrix(o, n.b2.x);
    matrix(o, n.w3.x);
    matrix(o, n.b3.x);
  }
  if (!o)
    throw std::runtime_error("write QFLOW failed");
}
} // namespace

int main(int argc, char **argv) {
  try {
    auto o = options(argc, argv);
    auto c = corpus(o.draws);
    std::map<int, int> mx;
    for (auto g : c.groups)
      mx[g.first] = std::max(mx[g.first], g.second);
    std::vector<int> tr, va;
    for (size_t i = 0; i < c.groups.size(); ++i)
      (c.groups[i].second == mx[c.groups[i].first] ? va : tr).push_back(i);
    auto rows = [&](const std::vector<int> &ix) {
      Matrix x(ix.size(), c.x.cols());
      for (int i = 0; i < x.rows(); ++i)
        x.row(i) = c.x.row(ix[i]);
      return x;
    };
    Matrix train = rows(tr), valid = rows(va);
    Flow f(c.x.cols(), o.layers, o.hidden, o.seed);
    f.mean = train.colwise().mean();
    f.scale = ((train.rowwise() - f.mean).array().square().colwise().sum() /
               float(train.rows() - 1))
                  .sqrt();
    Flow best = f;
    float bv = INFINITY, bt = INFINITY;
    int stale = 0, step = 0, done = 0;
    double gradient_check_error = 0.0;
    std::mt19937_64 rng(o.seed);
    for (int e = 1; e <= o.epochs; ++e) {
      std::shuffle(tr.begin(), tr.end(), rng);
      for (size_t s = 0; s < tr.size(); s += 128) {
        size_t z = std::min<size_t>(128, tr.size() - s);
        Matrix b(z, c.x.cols());
        for (size_t i = 0; i < z; ++i)
          b.row(i) = c.x.row(tr[s + i]);
        f.train_batch(b, ++step, 1e-3f);
        if (step == 1)
          gradient_check_error =
              check_gradients(f, b.topRows(std::min<int>(8, b.rows())));
      }
      float tn = f.forward(train), vn = f.forward(valid);
      done = e;
      if (vn < bv - 1e-4) {
        bv = vn;
        bt = tn;
        best = f;
        stale = 0;
      } else
        ++stale;
      if (e == 1 || e % 25 == 0)
        std::cout << "epoch=" << e << " training_nll=" << tn
                  << " validation_nll=" << vn << "\n";
      if (stale >= o.patience)
        break;
    }
    f = best;
    std::filesystem::create_directories(
        std::filesystem::path(o.output).parent_path());
    save(f, o.output);
    quadra::transport::DependencyFreeRealNVP check;
    check.load(o.output);
    double inverse_error = 0;
    std::normal_distribution<double> normal;
    for (int i = 0; i < 128; ++i) {
      std::vector<double> z(c.names.size());
      for (double &v : z)
        v = normal(rng);
      auto x = check.inverse(z), roundtrip = check.forward_noise(x);
      for (size_t j = 0; j < z.size(); ++j)
        inverse_error = std::max(inverse_error, std::abs(z[j] - roundtrip[j]));
    }
    if (!std::isfinite(inverse_error) || inverse_error > 1e-4)
      throw std::runtime_error("QFLOW round-trip validation failed");
    std::ofstream m(o.output + ".manifest");
    m << std::setprecision(17)
      << "manifest_version=1\nartifact_hash_fnv1a64=" << hash_file(o.output)
      << "\nassessment_fingerprint=" << fingerprint(o.fit)
      << "\ngeometry_fingerprint=" << fingerprint(o.geometry)
      << "\ndimension=" << c.names.size() << "\nparameter_names=";
    for (size_t i = 0; i < c.names.size(); ++i)
      m << (i ? "," : "") << c.names[i];
    m << "\narchitecture=qflow_realnvp_cpp;layers=" << o.layers
      << ";hidden=" << o.hidden
      << ";log_scale_limit=1.5\ntraining_draw_sources=";
    for (size_t i = 0; i < o.draws.size(); ++i)
      m << (i ? ";" : "") << o.draws[i];
    m << "\ntraining_draw_hashes_fnv1a64=";
    for (size_t i = 0; i < o.draws.size(); ++i)
      m << (i ? ";" : "") << hash_file(o.draws[i]);
    m << "\ntraining_draws=" << tr.size() << "\nvalidation_draws=" << va.size()
      << "\nvalidation_nll=" << bv << "\ntraining_nll=" << bt
      << "\ntraining_seed=" << o.seed << "\nepochs_completed=" << done
      << "\ninverse_max_error=" << inverse_error
      << "\ngradient_check_max_scaled_error=" << gradient_check_error << "\n";
    std::cout << "dependency-free RealNVP trained: training_nll=" << bt
              << " validation_nll=" << bv << " gap=" << bv - bt
              << " inverse_max_error=" << inverse_error
              << " gradient_error=" << gradient_check_error << "\n";
    return 0;
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
