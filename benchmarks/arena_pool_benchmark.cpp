#include "../core/memory/arena_pool.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

struct SmallNode {
    double value;
    double adjoint;
    int op;
    int left;
    int right;

    SmallNode(double value_, double adjoint_, int op_, int left_, int right_)
        : value(value_), adjoint(adjoint_), op(op_), left(left_), right(right_) {}
};

template <class F>
double time_ms(F&& f) {
    auto t0 = std::chrono::high_resolution_clock::now();
    f();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

int main() {
    const std::vector<std::size_t> Ns = {
        1000,
        10000,
        100000,
        1000000
    };

    const int repeats = 20;

    std::cout << "\nQuadra arena pool allocation benchmark\n\n";
    std::cout << std::setw(12) << "N"
              << std::setw(18) << "new/delete ms"
              << std::setw(18) << "arena ms"
              << std::setw(14) << "speedup"
              << std::setw(18) << "arena MB"
              << std::setw(12) << "blocks"
              << "\n";
    std::cout << std::string(92, '-') << "\n";

    for (std::size_t N : Ns) {
        double heap_ms = time_ms([&]() {
            for (int r = 0; r < repeats; ++r) {
                std::vector<SmallNode*> nodes;
                nodes.reserve(N);

                for (std::size_t i = 0; i < N; ++i) {
                    nodes.push_back(new SmallNode(static_cast<double>(i), 0.0, 1, -1, -1));
                }

                volatile double sink = 0.0;
                for (auto* node : nodes) {
                    sink += node->value;
                }

                for (auto* node : nodes) {
                    delete node;
                }
            }
        });

        quadra::ArenaPool arena(1u << 20);

        double arena_ms = time_ms([&]() {
            for (int r = 0; r < repeats; ++r) {
                arena.reset();

                std::vector<SmallNode*, quadra::ArenaAllocator<SmallNode*>>
                    nodes{quadra::ArenaAllocator<SmallNode*>(&arena)};
                nodes.reserve(N);

                for (std::size_t i = 0; i < N; ++i) {
                    nodes.push_back(arena.create<SmallNode>(static_cast<double>(i), 0.0, 1, -1, -1));
                }

                volatile double sink = 0.0;
                for (auto* node : nodes) {
                    sink += node->value;
                }
            }
        });

        const double speedup = arena_ms > 0.0 ? heap_ms / arena_ms : 0.0;
        const double arena_mb = static_cast<double>(arena.bytes_reserved()) / (1024.0 * 1024.0);

        std::cout << std::setw(12) << N
                  << std::setw(18) << std::fixed << std::setprecision(3) << heap_ms
                  << std::setw(18) << arena_ms
                  << std::setw(14) << speedup
                  << std::setw(18) << arena_mb
                  << std::setw(12) << arena.block_count()
                  << "\n";
    }

    std::cout << "\nInterpretation:\n";
    std::cout << "  If arena allocation is much faster here, allocator churn is worth pursuing.\n";
    std::cout << "  The next patch should wire ArenaPool into a had_quadra ADGraph wrapper,\n";
    std::cout << "  not globally rewrite had.h immediately.\n\n";

    return 0;
}
