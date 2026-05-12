#include "../core/autodiff/had_quadra_workspace.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

struct FakeADNode {
    double value;
    double adjoint;
    double hdot;
    int op;
    int left;
    int right;

    FakeADNode(double value_, int op_, int left_, int right_)
        : value(value_), adjoint(0.0), hdot(0.0), op(op_), left(left_), right(right_) {}
};

struct FakeHessianEdge {
    unsigned int index;
    double value;

    FakeHessianEdge(unsigned int index_, double value_)
        : index(index_), value(value_) {}
};

template <class F>
double time_ms(F&& f) {
    auto t0 = std::chrono::high_resolution_clock::now();
    f();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

static volatile double global_sink = 0.0;

void run_heap_eval(std::size_t nodes_per_eval, std::size_t edges_per_eval) {
    std::vector<FakeADNode*> nodes;
    nodes.reserve(nodes_per_eval);

    std::vector<FakeHessianEdge*> edges;
    edges.reserve(edges_per_eval);

    for (std::size_t i = 0; i < nodes_per_eval; ++i) {
        nodes.push_back(new FakeADNode(static_cast<double>(i), 1, -1, -1));
    }

    for (std::size_t i = 0; i < edges_per_eval; ++i) {
        edges.push_back(new FakeHessianEdge(static_cast<unsigned int>(i), 0.001 * static_cast<double>(i)));
    }

    double acc = 0.0;
    for (auto* node : nodes) acc += node->value + node->adjoint + node->hdot;
    for (auto* edge : edges) acc += edge->value;
    global_sink += acc;

    for (auto* edge : edges) delete edge;
    for (auto* node : nodes) delete node;
}

void run_workspace_eval(
    quadra::HadQuadraWorkspace& workspace,
    std::size_t nodes_per_eval,
    std::size_t edges_per_eval) {

    workspace.reset();

    auto nodes = workspace.make_vector_with_reserve<FakeADNode*>(nodes_per_eval);
    auto edges = workspace.make_vector_with_reserve<FakeHessianEdge*>(edges_per_eval);

    for (std::size_t i = 0; i < nodes_per_eval; ++i) {
        nodes.push_back(workspace.create<FakeADNode>(static_cast<double>(i), 1, -1, -1));
    }

    for (std::size_t i = 0; i < edges_per_eval; ++i) {
        edges.push_back(workspace.create<FakeHessianEdge>(
            static_cast<unsigned int>(i),
            0.001 * static_cast<double>(i)));
    }

    double acc = 0.0;
    for (auto* node : nodes) acc += node->value + node->adjoint + node->hdot;
    for (auto* edge : edges) acc += edge->value;
    global_sink += acc;
}

int main() {
    struct Case {
        std::size_t G;
        std::size_t m;
        std::size_t evals;
    };

    const std::vector<Case> cases = {
        {25, 20, 50},
        {50, 20, 50},
        {100, 20, 50},
        {250, 20, 50},
        {500, 20, 50},
        {1000, 20, 50}
    };

    std::cout << "\nHadQuadraWorkspace AD-lifecycle allocation benchmark\n\n";
    std::cout << std::setw(8) << "G"
              << std::setw(8) << "m"
              << std::setw(10) << "evals"
              << std::setw(16) << "nodes/eval"
              << std::setw(16) << "edges/eval"
              << std::setw(16) << "heap ms"
              << std::setw(18) << "workspace ms"
              << std::setw(14) << "speedup"
              << std::setw(16) << "reserved MB"
              << std::setw(10) << "blocks"
              << "\n";
    std::cout << std::string(130, '-') << "\n";

    for (const auto& c : cases) {
        const std::size_t nodes_per_eval = c.G * (c.m + 8);
        const std::size_t edges_per_eval = 3 * c.G + 2 * c.m;

        double heap_ms = time_ms([&]() {
            for (std::size_t e = 0; e < c.evals; ++e) {
                run_heap_eval(nodes_per_eval, edges_per_eval);
            }
        });

        quadra::HadQuadraWorkspace workspace(1u << 20);

        double workspace_ms = time_ms([&]() {
            for (std::size_t e = 0; e < c.evals; ++e) {
                run_workspace_eval(workspace, nodes_per_eval, edges_per_eval);
            }
        });

        const double speedup = workspace_ms > 0.0 ? heap_ms / workspace_ms : 0.0;
        const double reserved_mb =
            static_cast<double>(workspace.bytes_reserved()) / (1024.0 * 1024.0);

        std::cout << std::setw(8) << c.G
                  << std::setw(8) << c.m
                  << std::setw(10) << c.evals
                  << std::setw(16) << nodes_per_eval
                  << std::setw(16) << edges_per_eval
                  << std::setw(16) << std::fixed << std::setprecision(3) << heap_ms
                  << std::setw(18) << workspace_ms
                  << std::setw(14) << speedup
                  << std::setw(16) << reserved_mb
                  << std::setw(10) << workspace.block_count()
                  << "\n";
    }

    std::cout << "\nResult checksum sink: " << global_sink << "\n\n";

    std::cout << "Next integration target if this is strong:\n";
    std::cout << "  1. Add optional HadQuadraWorkspace* argument to evaluate_random_effect_hessian().\n";
    std::cout << "  2. Move temporary std::vector work buffers into workspace.make_vector_with_reserve<T>().\n";
    std::cout << "  3. Reset workspace once per evaluation, not per inner loop.\n";
    std::cout << "  4. Only after that, consider a narrow had::NewAReal arena-backed path.\n\n";

    return 0;
}
