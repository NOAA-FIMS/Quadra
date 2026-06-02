#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

// Isolated prototype for forward replay metadata.
//
// This test does NOT modify production had_quadra.hpp yet.
// It demonstrates the minimum metadata needed for reusable forward replay.

namespace prototype {

enum class OpCode {
  Constant,
  Independent,
  Add,
  Subtract,
  Multiply,
  Divide,
  Negate
};

struct Node {
  OpCode op_m = OpCode::Constant;
  int left_m = -1;
  int right_m = -1;
  double value_m = 0.0;
};

class ReplayGraph {
public:
  int constant(double value) {
    nodes_m.push_back({OpCode::Constant, -1, -1, value});
    return static_cast<int>(nodes_m.size()) - 1;
  }

  int independent(double value) {
    nodes_m.push_back({OpCode::Independent, -1, -1, value});
    return static_cast<int>(nodes_m.size()) - 1;
  }

  int add(int left, int right) {
    nodes_m.push_back({OpCode::Add, left, right, value(left) + value(right)});
    return static_cast<int>(nodes_m.size()) - 1;
  }

  int subtract(int left, int right) {
    nodes_m.push_back(
        {OpCode::Subtract, left, right, value(left) - value(right)});
    return static_cast<int>(nodes_m.size()) - 1;
  }

  int multiply(int left, int right) {
    nodes_m.push_back(
        {OpCode::Multiply, left, right, value(left) * value(right)});
    return static_cast<int>(nodes_m.size()) - 1;
  }

  int divide(int left, int right) {
    nodes_m.push_back(
        {OpCode::Divide, left, right, value(left) / value(right)});
    return static_cast<int>(nodes_m.size()) - 1;
  }

  int negate(int left) {
    nodes_m.push_back({OpCode::Negate, left, -1, -value(left)});
    return static_cast<int>(nodes_m.size()) - 1;
  }

  void set_value(int id, double value) {
    nodes_m.at(static_cast<size_t>(id)).value_m = value;
  }

  double value(int id) const {
    return nodes_m.at(static_cast<size_t>(id)).value_m;
  }

  void forward() {
    for (size_t i = 0; i < nodes_m.size(); ++i) {
      Node &node = nodes_m[i];

      switch (node.op_m) {
      case OpCode::Constant:
      case OpCode::Independent:
        break;

      case OpCode::Add:
        node.value_m = value(node.left_m) + value(node.right_m);
        break;

      case OpCode::Subtract:
        node.value_m = value(node.left_m) - value(node.right_m);
        break;

      case OpCode::Multiply:
        node.value_m = value(node.left_m) * value(node.right_m);
        break;

      case OpCode::Divide:
        node.value_m = value(node.left_m) / value(node.right_m);
        break;

      case OpCode::Negate:
        node.value_m = -value(node.left_m);
        break;

      default:
        throw std::runtime_error("Unknown replay opcode.");
      }
    }
  }

  size_t size() const { return nodes_m.size(); }

private:
  std::vector<Node> nodes_m;
};

} // namespace prototype

int main() {
  std::cout << "forward replay design prototype\n";

  prototype::ReplayGraph graph;

  int x = graph.independent(2.0);
  int one = graph.constant(1.0);
  int x2 = graph.multiply(x, x);
  int y = graph.add(x2, one);

  std::cout << "initial y = " << graph.value(y) << "\n";

  if (std::abs(graph.value(y) - 5.0) > 1e-12) {
    std::cerr << "FAIL: initial replay value mismatch\n";
    return 1;
  }

  graph.set_value(x, 3.0);
  graph.forward();

  std::cout << "mutated y after forward replay = " << graph.value(y) << "\n";

  if (std::abs(graph.value(y) - 10.0) > 1e-12) {
    std::cerr << "FAIL: replay value mismatch\n";
    return 1;
  }

  graph.set_value(x, -4.0);
  graph.forward();

  std::cout << "second mutated y after forward replay = " << graph.value(y)
            << "\n";

  if (std::abs(graph.value(y) - 17.0) > 1e-12) {
    std::cerr << "FAIL: second replay value mismatch\n";
    return 1;
  }

  std::cout << "node count = " << graph.size() << "\n";
  std::cout << "PASS\n";

  return 0;
}
