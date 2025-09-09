cmake_minimum_required(VERSION 3.28.3)
project(untitled)

set(CMAKE_CXX_STANDARD 20)

find_package(benchmark QUIET)

if(NOT benchmark_FOUND)
    include(FetchContent)

    FetchContent_Declare(
            googlebenchmark
            GIT_REPOSITORY https://github.com/google/benchmark.git
            GIT_TAG v1.8.3
    )

    set(BENCHMARK_ENABLE_TESTING OFF)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF)

    FetchContent_MakeAvailable(googlebenchmark)
endif ()

option(MAIN "Build main exec" ON)
option(BENCH "Build benchmark" ON)

if (MAIN)
    add_executable(main main.cpp solver.cpp)
endif ()

if (BENCH)
    add_executable(bench bench.cpp solver.cpp solver.h)
    target_link_libraries(bench PRIVATE benchmark::benchmark)
endif ()