#include "matrices/tridiagonal_matrix.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::matrices);

void TridiagonalMatrix::set(int i, int j, double value)
{
    assert(i >= 0);
    assert(j >= 0);
    assert(i < size);
    assert(j < size);
    assert (!std::isnan(value));

    if (i == j)
    {
        b[i] = value; // Главная диагональ
        return;
    }
    if (i == j + 1)
    {
        a[j] = value; // Нижняя диагональ
        return;
    }
    if (i == j - 1)
    {
        c[i] = value; // Верхняя диагональ
        return;
    }
    throw std::runtime_error("Only tridiagonal elements can be non-zero");
}

std::vector<double> TridiagonalMatrix::solve(std::vector<double> const &d) const
{
    if (d.size() != static_cast<size_t>(size))
    {
        throw std::invalid_argument("Vector size must match matrix size");
    }

    std::vector<double> x(size);
    std::vector<double> alpha(size - 1);
    std::vector<double> beta(size);

    // Прямой ход прогонки
    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];

    for (int i = 1; i < size - 1; ++i)
    {
        double denominator = b[i] + a[i - 1] * alpha[i - 1];
        alpha[i] = -c[i] / denominator;
        beta[i] = (d[i] - a[i - 1] * beta[i - 1]) / denominator;
    }

    // Обратный ход прогонки
    beta[size - 1] = (d[size - 1] - a[size - 2] * beta[size - 2]) / (b[size - 1] + a[size - 2] * alpha[size - 2]);
    x[size - 1] = beta[size - 1];

    for (int i = size - 2; i >= 0; --i)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    return x;
}

PLSAPR_END_NAMESPACE(); // plastinka_sapr::matrices
