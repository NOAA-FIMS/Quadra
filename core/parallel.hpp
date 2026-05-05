#ifndef PARALLEL_HPP
#define PARALLEL_HPP
#pragma once

#include <vector>
#include <thread>
#include <algorithm>

namespace quadra
{

    // Generic parallel reduce over index range
    template <typename Func>
    double parallel_reduce(
        int n_items,
        Func &&func)
    {
        unsigned int n_threads = std::thread::hardware_concurrency();
        if (n_threads == 0)
            n_threads = 4;

        int chunk = (n_items + n_threads - 1) / n_threads;

        std::vector<std::thread> threads;
        std::vector<double> results(n_threads, 0.0);

        for (unsigned int t = 0; t < n_threads; ++t)
        {

            int start = t * chunk;
            int end = std::min(start + chunk, n_items);

            if (start >= end)
                break;

            threads.emplace_back([&, t, start, end]()
                                 {

            double local = 0.0;

            for (int i = start; i < end; ++i)
                local += func(i);

            results[t] = local; });
        }

        for (auto &th : threads)
            th.join();

        double total = 0.0;
        for (double v : results)
            total += v;

        return total;
    }

} // namespace pelagia

#endif // PARALLEL_HPP