#include "matrices/eigen_sparce_matrix.hpp"

#include <cassert>
#include <cmath>

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::matrices);

void BaseSparceMatrix::set(int i, int j, double value)
{
    assert(i >= 0);
    assert(j >= 0);
    assert(i < size);
    assert(j < size);
    assert(!std::isnan(value));

    matrix.insert(i, j) = value;
}

std::vector<double> SparceMatrixLU::solve(const std::vector<double> &x) const
{
    auto matrixCopy = matrix;
    Eigen::VectorXd vector = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; i++)
    {
        vector(i) = x[i];
    }

    matrixCopy.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(matrix);
    solver.factorize(matrix);

    const auto solution = solver.solve(vector);
    auto result = std::vector<double>(size);
    for (int i = 0; i < size; i++)
    {
        result[i] = solution(i);
    }
    return result;
}

PLSAPR_END_NAMESPACE(); // plastinka_sapr::matrices
