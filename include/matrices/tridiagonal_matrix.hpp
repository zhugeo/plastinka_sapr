#pragma once

#include <vector>

#include "base_matrix.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::matrices);

class TridiagonalMatrix : public BaseMatrix
{
private:
    std::vector<double> a; // Нижняя диагональ
    std::vector<double> b; // Главная диагональ
    std::vector<double> c; // Верхняя диагональ

public:
    TridiagonalMatrix(int size) : BaseMatrix(size),
                                  a(std::vector<double>(size - 1, 0.0)),
                                  b(std::vector<double>(size, 0.0)),
                                  c(std::vector<double>(size - 1, 0.0)) {}
    virtual void set(int i, int j, double value) override;
    virtual std::vector<double> solve(const std::vector<double> &x) const override;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr::matrices
