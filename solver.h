//
// Created by ih284669 on 05/09/2025.
//

#ifndef UNTITLED_SOLVER_H
#define UNTITLED_SOLVER_H
#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numeric>

struct CSR {
    std::vector<size_t> cols;
    std::vector<double> values;
    std::vector<size_t> rowsPtr;
};

class Solver {
    int gridSize;
    double dt;
    double dx;
    double k;
    double D;
    double Tref;
    size_t nb_iterations = 0;
    std::vector<double> oldGrid;
    std::vector<double> newGrid;

public:
    Solver(int gridSize, double dt, double dx, double diffusion, const std::vector<double>& initGrid, double k, double T);

    void solveNexplicit(const size_t n);

    std::vector<double> multiplyScalar(std::vector<double>& result, const double a);

    std::vector<double> matrixMultiply(const std::vector<double>& mat, const std::vector<double>& v);

    std::vector<double> matrixMultiplySparsev1(const CSR& mat, const std::vector<double>& v);
    std::vector<double> matrixMultiplySparsev2(const CSR& mat, const std::vector<double>& v);

    void fillA(std::vector<double>& A);

    void fillAsparse(CSR& A);

    void solveImplicit(size_t max_iterGC, size_t n, bool sparse);

    static double refSolution(const double k, const double dt, const double dx, const double i, const double n);

    std::vector<double> getGrid();

    std::vector<double> getRef();

    double getL2Error();
};
#endif //UNTITLED_SOLVER_H