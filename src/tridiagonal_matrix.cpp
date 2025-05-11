#include "tridiagonal_matrix.hpp"

void TridiagonalMatrix::set(int i, int j, double value)
{

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
    if (d.size() != static_cast<size_t>(n))
    {
        throw std::invalid_argument("Vector size must match matrix size");
    }

    std::vector<double> x(n);
    std::vector<double> alpha(n - 1);
    std::vector<double> beta(n);

    // Прямой ход прогонки
    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];

    for (int i = 1; i < n - 1; ++i)
    {
        double denominator = b[i] + a[i - 1] * alpha[i - 1];
        alpha[i] = -c[i] / denominator;
        beta[i] = (d[i] - a[i - 1] * beta[i - 1]) / denominator;
    }

    // Обратный ход прогонки
    beta[n - 1] = (d[n - 1] - a[n - 2] * beta[n - 2]) / (b[n - 1] + a[n - 2] * alpha[n - 2]);
    x[n - 1] = beta[n - 1];

    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    return x;
}
