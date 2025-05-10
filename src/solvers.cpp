#include "solvers.hpp"

#include <cmath>

#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

namespace
{
    bool isConsistent(const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
    {
        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
        int rankA = lu_decomp.rank();

        Eigen::MatrixXd Ab(A.rows(), A.cols() + 1);
        Ab << A, b;

        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_Ab(Ab);
        int rankAb = lu_decomp_Ab.rank();

        return (rankA == rankAb);
    }
} // namespace

Solution solveImplicit(const Model &model, const std::shared_ptr<Grid> grid, double dt, int numIters, double initT)
{
    auto T = std::vector<std::vector<double>>(numIters);
    int N = grid->innerNodes.size() + grid->outerNodes.size();

    // Расчёт начального состояния
    T[0] = std::vector<double>(N, initT);
    for (int i = 0; i < grid->outerNodes.size(); i++)
    {
        auto &node = grid->outerNodes[i];
        auto idx = i + grid->innerNodes.size();
    }

    // Расчёт последующих итераций
    for (int iter = 1; iter < numIters; iter++)
    {
        Eigen::SparseMatrix<double> matrix(N, N);
        Eigen::VectorXd vector = Eigen::VectorXd::Zero(N);
        for (int nodeIndex = 0; nodeIndex < grid->innerNodes.size(); nodeIndex++)
        {
            const auto node = grid->innerNodes[nodeIndex];

            const auto leftNode = node->left.lock();
            const auto rightNode = node->right.lock();
            const auto topNode = node->top.lock();
            const auto bottomNode = node->bottom.lock();

            const int leftNodeIndex = grid->getNodeIndex(node->left.lock());
            const int rightNodeIndex = grid->getNodeIndex(node->right.lock());
            const int topNodeIndex = grid->getNodeIndex(node->top.lock());
            const int bottomNodeIndex = grid->getNodeIndex(node->bottom.lock());

            double mu_x = rightNode->muValue;
            if (mu_x == 1)
            {
                mu_x = 1 / leftNode->muValue;
            }
            double mu_y = topNode->muValue;
            if (mu_y == 1)
            {
                mu_y = 1 / bottomNode->muValue;
            }

            double K = 2 * model.lambda * dt / model.c / model.rho;
            double Kx = 1 / (mu_x * (mu_x + 1) * grid->xStep * grid->xStep);
            double Ky = 1 / (mu_y * (mu_y + 1) * grid->yStep * grid->yStep);

            matrix.insert(nodeIndex, nodeIndex) = (-K * Kx * (mu_x + 1) - K * Ky * (mu_y + 1) - 1);
            assert(!std::isnan(matrix.coeff(nodeIndex, nodeIndex)));
            matrix.insert(nodeIndex, leftNodeIndex) = K * Kx * mu_x;
            matrix.insert(nodeIndex, rightNodeIndex) = K * Kx;
            matrix.insert(nodeIndex, topNodeIndex) = K * Ky;
            matrix.insert(nodeIndex, bottomNodeIndex) = K * Ky * mu_y;
            vector(nodeIndex) = -T[iter - 1][nodeIndex];
        }
        for (int i = 0; i < grid->outerNodes.size(); i++)
        {
            const auto node = grid->outerNodes[i];
            const auto parentNode = node->parent;
            const auto nodeIndex = grid->getNodeIndex(node);
            const auto parentNodeIndex = grid->getNodeIndex(node->parent.lock());

            const auto borderType = node->type;
            const auto relevantStep = ((node->side == OuterNodeSide::Left || node->side == OuterNodeSide::Right) ? grid->xStep : grid->yStep);
            assert(node->muValue > 0);
            assert(node->muValue <= 1);
            assert(nodeIndex != parentNodeIndex);
            if (borderType == BorderType::constFlow)
            {
                matrix.insert(nodeIndex, parentNodeIndex) = -1;
                matrix.insert(nodeIndex, nodeIndex) = 1;
                vector(nodeIndex) = node->borderValue * node->muValue * relevantStep;
            }
            else if (borderType == BorderType::constTemperature)
            {
                matrix.insert(nodeIndex, nodeIndex) = 1;
                vector(nodeIndex) = node->borderValue;
            }
            else if (borderType == BorderType::convection)
            {
                matrix.insert(nodeIndex, nodeIndex) = 1 - node->borderValue * node->muValue * relevantStep;
                matrix.insert(nodeIndex, parentNodeIndex) = -1;
                vector(nodeIndex) = 0;
            }
        }

        // В целях дебага можно перейти на плотные матрицы и попробовать эти проверки
        // const auto matrixDeterminant = matrix.determinant();
        // assert(matrixDeterminant != 0);
        // assert(!std::isnan(matrixDeterminant));
        // assert(isConsistent(matrix, vector));

        matrix.makeCompressed();

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(matrix);
        solver.factorize(matrix);

        const auto solution = solver.solve(vector);

        T[iter] = std::vector<double>(N);
        for (int i = 0; i < N; i++)
        {
            T[iter][i] = solution(i);
        }
    }

    Solution sol;
    sol.T = T;
    sol.timeStep = dt;
    sol.grid = grid;
    return sol;
}
