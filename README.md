# MiniHeatSolver

## Getting Started

To use this solver see [Installation](#installation) and [Usage](#usage).
From here, you will be able to use main file to play around with the parameter of the solver or
to run benchmark to compare Explicit and Implicit resolution.

## Installation

This section describe how you can install MiniHeatSolver with git and how you build and compile it with cmake.
You'll need git and cmake installed.

```bash
# Clone the repo
$ git clone https://github.com/Input-Fr/MiniHeatSolver.git
# Go to the root directory
$ cd MiniHeatSolver
# Build both the main and the benchmark
$ cmake -DMAIN=ON -DBENCH=ON -B build
# Go to the build
$ cd build
# Compile the whole project
$ make
```

This compile the `main` and `bench` binaries.
The benchmark uses google benchmark, if you have already installed it, nothing more will be installed,
otherwise it will be temporarily installed in the build. 

> [!TIP]
> You can choose to not build the bench or the main (to build only one of the two binaries) by setting `-DBENCH=OFF` or `-DMAIN=OFF` while generating the build with `cmake`.


## Usage

A basic benchmark have been built to compare explicit and implicit methods.
A main is also available to use this mini solver freely.

```bash
# Execute base benchmark
$ ./bench 
Solving 1D Heat equations     CPU           L2 Error         MaxIter_CG         Time Step         Iterations         Time per Iteration
---------------------------------------------------------------------------------------------------------------------------------------
Explicit_Dense/real_time     1327.577651ms  9.07598e-11      0                   4.5e-07              222222               0.00597411
Implicit_Sparse/real_time    894.200393ms   5.79853e-05      100                 0.000225             444                  2.01396
# Execute base main
$ ./main 
L2 error : 5.99904e-32
Computed grid : [0,0.587785,0.951057,0.951057,0.587785,1.22465e-16,-0.587785,-0.951057,-0.951057,-0.587785,0]
Reference : [0,0.587785,0.951057,0.951057,0.587785,1.22465e-16,-0.587785,-0.951057,-0.951057,-0.587785,-2.44929e-16]
```
