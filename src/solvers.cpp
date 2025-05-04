#include "solvers.hpp"
#include "Eigen/Sparse"
#include "Eigen/Dense"

Solution solveImplicit(const Model &model, const Grid &grid, double dt, int numIters, double initT)
{
    auto T = std::vector<std::vector<double>>();
    int N = grid.innerNodes.size() + grid.outerNodes.size();

    // Расчёт начального состояния
    T.push_back(std::vector<double>(N, initT));
    for (int i = 0; i < grid.outerNodes.size(); i++)
    {
        auto &node = grid.outerNodes[i];
        auto idx = i + grid.innerNodes.size();
    }

    // Расчёт последующих итераций
    for (int iter = 1; iter < numIters; iter++)
    {
        Eigen::MatrixXd matrix(N, N);
        Eigen::VectorXd vector(N);
        for (int i = 0; i < grid.innerNodes.size(); i++)
        {
            auto &node = grid.innerNodes[i];
            matrix(i, i) = -static_cast<double>(2) / std::pow(grid.xStep, 2) - static_cast<double>(2) / std::pow(grid.yStep, 2) + static_cast<double>(2) / dt;
            matrix(i, grid.getNodeIndex(node->left.lock())) = std::pow(grid.xStep, -1);
            matrix(i, grid.getNodeIndex(node->right.lock())) = std::pow(grid.xStep, -1);
            matrix(i, grid.getNodeIndex(node->top.lock())) = std::pow(grid.yStep, -1);
            matrix(i, grid.getNodeIndex(node->bottom.lock())) = std::pow(grid.yStep, -1);
            vector(i) = T.back()[i] / dt;
        }
        for (int i = 0; i < grid.outerNodes.size(); i++)
        {
            auto &node = grid.outerNodes[i];
            auto row = i + grid.innerNodes.size();
        }
    }

    Solution sol;
    sol.T = T;
    return sol;
}
