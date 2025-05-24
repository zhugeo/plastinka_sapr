#pragma once

#include <vector>

#include <Eigen/SparseLU>

#include "base_matrix.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::matrices);

class BaseSparceMatrix : public BaseMatrix
{
protected:
    Eigen::SparseMatrix<double> matrix;

public:
    virtual void set(int i, int j, double value) final override;
    virtual std::vector<double> solve(const std::vector<double> &x) const = 0;

    BaseSparceMatrix(int size) : BaseMatrix(size), matrix(size, size) {}
};

class SparceMatrixLU final : public BaseSparceMatrix
{
public:
    virtual std::vector<double> solve(const std::vector<double> &x) const;

    SparceMatrixLU(int size) : BaseSparceMatrix(size) {}
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr::matrices
