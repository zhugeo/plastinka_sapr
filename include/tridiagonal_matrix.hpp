#pragma once

#include <vector>
#include <stdexcept>

#include "base.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr);

class TridiagonalMatrix
{
private:
    int n;
    std::vector<double> a; // Нижняя диагональ
    std::vector<double> b; // Главная диагональ
    std::vector<double> c; // Верхняя диагональ

public:
    TridiagonalMatrix(int size)
        : n(size),
          a(std::vector<double>(n - 1, 0.0)),
          b(std::vector<double>(n, 0.0)),
          c(std::vector<double>(n - 1, 0.0))
    {
    }
    void set(int i, int j, double value);
    std::vector<double> solve(const std::vector<double> &d) const;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr
