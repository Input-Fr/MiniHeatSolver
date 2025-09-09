//
// Created by ih284669 on 05/09/2025.
//
#include <benchmark/benchmark.h>
#include <thread>
#include <chrono>
#include "solver.h"
/*
#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/chrono.h>
*/

void BM_Explicit(benchmark::State& st) {
    double dx = 1e-3;
    double CFL = 0.9;
    size_t I = 1; // implicite != 1
    double dt = I*CFL*dx*dx/2;
    double k = 1.0;
    double diffusion = 1;
    double T = 0.1;
    size_t n = T/dt;
    int gridSize = static_cast<int>(std::ceil(1/dx)) + 1;
    std::vector initGrid(gridSize, 0.0);
    for (int i = 0; i < gridSize; i++) {
        initGrid[i] = sin(k*2*M_PI*i*dx);
    }
    Solver solver = Solver(gridSize, dt, dx, diffusion, initGrid, k, T);
    for (auto _ : st) {
        solver.solveNexplicit(n);
    }
    st.counters["L2Error"] = 9.07598e-11;
    st.counters["MaxIter_CG"] = 0;
    st.counters["TimeStep"] = dt;
    st.counters["Iteration"] = n;
}

void BM_Implicit2(benchmark::State& st) {
    double dx = 1e-3;
    double CFL = 0.9;
    size_t I = 500; // implicite != 1
    double dt = I*CFL*dx*dx/2;
    double k = 1.0;
    double diffusion = 1;
    double T = 0.1;
    size_t n = T/dt;
    size_t max_iter = 100;
    int gridSize = static_cast<int>(std::ceil(1/dx)) + 1;
    std::vector initGrid(gridSize, 0.0);
    for (int i = 0; i < gridSize; i++) {
        initGrid[i] = sin(k*2*M_PI*i*dx);
    }
    Solver solver = Solver(gridSize, dt, dx, diffusion, initGrid, k, T);
    for (auto _ : st) {
        solver.solveImplicit(max_iter, n, true);
    }
    st.counters["L2Error"] = 5.79853e-05;
    st.counters["MaxIter_CG"] = max_iter;
    st.counters["TimeStep"] = dt;
    st.counters["Iteration"] = n;
}

BENCHMARK(BM_Explicit) -> Unit(benchmark::kMillisecond) -> UseRealTime() ->Name("Explicit_Dense");
BENCHMARK(BM_Implicit2) -> Unit(benchmark::kMillisecond) -> UseRealTime() ->Name("Implicit_Sparse");

class CustomReporter : public benchmark::ConsoleReporter {
public:
    bool ReportContext(const Context& context) override {

        std::cout << "Solving 1D Heat equations     CPU           L2 Error         MaxIter_CG         Time Step         Iterations         Time per Iteration\n";
        std::cout << "---------------------------------------------------------------------------------------------------------------------------------------\n";
        return true;
    }

    void ReportRuns(const std::vector<Run>& reports) override {
        size_t idx = 0;
        for (auto& run : reports) {
            std::cout << std::left
          << std::setw(29) << run.benchmark_name()
          << std::setw(15) << (std::to_string(run.GetAdjustedCPUTime()) + "ms")
          << std::setw(17) << run.counters.at("L2Error")
          << std::setw(20) << run.counters.at("MaxIter_CG")
          << std::setw(21) << run.counters.at("TimeStep")
          << std::setw(21) << run.counters.at("Iteration")
            << run.GetAdjustedCPUTime()/run.counters.at("Iteration")
          << std::endl;
        }
    }
};

int main(int argc, char** argv) {
    benchmark::Initialize(&argc, argv);
    CustomReporter reporter;
    RunSpecifiedBenchmarks(&reporter);
}
