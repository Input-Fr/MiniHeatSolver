//
// Created by ih284669 on 05/09/2025.
//
#include "solver.h"

Solver::Solver(int gridSize, double dt, double dx, double diffusion, const std::vector<double>& initGrid, double k, double T)
    : D(diffusion), gridSize(gridSize), dt(dt), dx(dx), oldGrid(initGrid), newGrid(initGrid), k(k), Tref(T) {}

void Solver::solveNexplicit(const size_t n) {
    double coef = D*(dt/(dx*dx));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 1; j < gridSize-1; j++) {
            newGrid[j] = oldGrid[j] + coef*(oldGrid[(j-1)] - 2*oldGrid[j] + oldGrid[(j+1)]);
        }
        oldGrid = newGrid;
    }
    nb_iterations += n;
}

std::vector<double> Solver::multiplyScalar(std::vector<double>& result, const double a) {
    std::vector temp(result);
    for (size_t i = 0; i < result.size(); i++) {
        temp[i] = result[i]*a;
    }
    return temp;
}

std::vector<double> Solver::matrixMultiply(const std::vector<double>& mat, const std::vector<double>& v) {
    std::vector rslt(gridSize, 0.0);
    for (int y = 0; y < gridSize; y++) {
        for (int x = 0; x < gridSize; x++) {
            rslt[y] += mat[x + y * gridSize] * v[x];
        }
    }
    return rslt;
}

std::vector<double> Solver::matrixMultiplySparsev1(const CSR& mat, const std::vector<double>& v) {
    std::vector rslt(gridSize, 0.0);
    double tmp = 0.0;
    size_t row = 0;
    for (int y = 0; y < mat.values.size(); y++) {
        if (y+1 < mat.cols.size() && mat.cols[y] < mat.cols[y+1]) {
            tmp+=mat.values[y]*v[mat.cols[y]];
        }
        else {
            tmp+=mat.values[y]*v[mat.cols[y]];
            rslt[row] = tmp;
            tmp = 0.0;
            row++;
        }
    }
    return rslt;
}

std::vector<double> Solver::matrixMultiplySparsev2(const CSR& mat, const std::vector<double>& v) {
    std::vector rslt(gridSize, 0.0);
    double tmp = 0.0;
    size_t row = 0;
    for (int y = 0; y < mat.cols.size(); y++) {
        if (y+1 != mat.rowsPtr[row+1]) {
            tmp+=mat.values[y]*v[mat.cols[y]];
        }
        else {
            tmp+=mat.values[y]*v[mat.cols[y]];
            rslt[row] = tmp;
            tmp = 0.0;
            row++;
        }
    }
    return rslt;
}

void Solver::fillA(std::vector<double>& A) {
    double coef = D*dt/(dx*dx);
    double diag = 1.0+2.0*coef;
    A[0] = 1;
    for (size_t y = 1; y < gridSize-1; y++) {
        A[y-1 + y*gridSize] = -coef;
        A[y + y*gridSize] = diag;
        A[y+1 + y*gridSize] = -coef;
    }
    A[gridSize*gridSize-1] = 1;
}

void Solver::fillAsparse(CSR& A) {
    double coef = D*dt/(dx*dx);
    double diag = 1.0+2.0*coef;
    A.values = std::vector<double>(0);
    A.cols = std::vector<size_t>(0);
    A.rowsPtr = std::vector<size_t>(0);

    A.cols.emplace_back(0);
    A.rowsPtr.emplace_back(0);
    A.values.emplace_back(1);

    for (int i = 1; i < gridSize-1; i++) {
        A.rowsPtr.emplace_back(A.values.size());
        A.cols.emplace_back(i-1);
        A.cols.emplace_back(i);
        A.cols.emplace_back(i+1);
        A.values.emplace_back(-coef);
        A.values.emplace_back(diag);
        A.values.emplace_back(-coef);
    }
    A.rowsPtr.emplace_back(A.values.size());
    A.cols.emplace_back(gridSize-1);
    A.values.emplace_back(1);

    //A.rowsPtr.emplace_back(A.cols.size()-1);
    // no need for this particular case
}


void Solver::solveImplicit(size_t max_iterGC, size_t n, bool sparse) {
    oldGrid[0] = 0.0;
    oldGrid[gridSize-1] = 0.0;
    newGrid = oldGrid;
    std::vector A(gridSize*gridSize, 0.0);
    CSR A_csr;
    sparse ? fillAsparse(A_csr) : fillA(A);
    double T = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (T+dt > Tref) {
            double saveDt = dt;
            dt = Tref-T;
            T += dt;
            sparse ? fillAsparse(A_csr) : fillA(A);
            dt = saveDt;
        }
        else
            T += dt;
        std::vector<double> r(oldGrid.size());
        std::vector Ax0 = sparse ? matrixMultiplySparsev2(A_csr,oldGrid) : matrixMultiply(A, oldGrid);
        std::transform(oldGrid.begin(), oldGrid.end(), Ax0.begin(), r.begin(), std::minus());
        std::vector p(r); // p0 = r0
        for (size_t idx = 0; idx < max_iterGC; idx++) {
            std::vector<double> Ap = sparse ? matrixMultiplySparsev2(A_csr,p) : matrixMultiply(A, p);
            double p_Ap = std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0);
            double rTr = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
            double alpha = rTr/p_Ap;
            std::vector<double> alP = multiplyScalar(p,alpha);
            std::transform(oldGrid.begin(), oldGrid.end(), alP.begin(), oldGrid.begin(), std::plus());
            std::transform(r.begin(), r.end(), multiplyScalar(Ap,alpha).begin(), r.begin(), std::minus());

            if (*std::ranges::max_element(r) < 1e-12) {
                break;
            }
            double beta = std::inner_product(r.begin(), r.end(), r.begin(), 0.0)/rTr;
            std::transform(r.begin(), r.end(), multiplyScalar(p,beta).begin(), p.begin(), std::plus());
        }
        oldGrid[0] = 0;
        oldGrid[gridSize-1] = 0;
    }
    newGrid = oldGrid;
    nb_iterations = n;
}

double Solver::refSolution(const double k, const double dt, const double dx, const double i, const double n) {
    return std::exp(-(k*2*M_PI)*(k*2*M_PI)*n*dt)*sin(k*2*M_PI*dx*i);
}

std::vector<double> Solver::getGrid() {
    return newGrid;
}

std::vector<double> Solver::getRef() {
    std::vector<double> ref(gridSize);
    for (size_t i = 0; i < gridSize; i++) {
        ref[i] = refSolution(k, dt, dx, static_cast<double>(i), static_cast<double>(nb_iterations));
    }
    return ref;
}

double L2(std::vector<double>& ref, std::vector<double>& pred) {
    double l2Error = 0.0;
    for (size_t i = 0; i < ref.size(); i++) {
        l2Error += std::pow(ref[i] - pred[i], 2);
    }
    return l2Error;
}

double Solver::getL2Error() {
    std::vector<double> ref = getRef();
    std::vector<double> pred = getGrid();
    return L2(ref, pred);
}