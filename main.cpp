#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numeric>
#include "solver.h"

void printVector(const std::vector<double>& v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << v[i] << (i == v.size() - 1 ? "" : ",");
    }
}

int main() {
    double dx = 1e-1;
    double CFL = 0.9;
    size_t I = 800; // implicite != 1
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
    solver.solveImplicit(100, n, true);


    // PRINT L2 ERROR
    std::cout << "L2 error : " << solver.getL2Error() << std::endl;

    // PRINT COMPUTED GRID
    std::vector<double> g = solver.getGrid();
    std::cout << "Computed grid : [";
    printVector(g);
    std::cout << ']' << std::endl;

    // PRINT REF GRID
    std::vector<double> r = solver.getRef();
    std::cout << "Reference : [";
    printVector(r);
    std::cout << ']' << std::endl;

    return 0;
}