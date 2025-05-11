#include "solvers.hpp"

#include <cmath>
#include <set>
#include <iostream>

#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include "tridiagonal_matrix.hpp"

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
    int N = grid->getNodeCount();

    // Расчёт начального состояния
    T[0] = std::vector<double>(N, initT);

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

Solution solveImplicitFast(const Model &model, const std::shared_ptr<Grid> grid, double dt, int numIters, double initT)
{
    auto T = std::vector<std::vector<double>>(numIters);
    int N = grid->getNodeCount();

    // Расчёт начального состояния
    T[0] = std::vector<double>(N, initT);

    // Расчёт последующих итераций
    for (int iter = 1; iter < numIters; iter++)
    {
        T[iter] = std::vector<double>(N);
        /** Горизонтальный проход */
        std::set<int> unindexedOuterNodes;
        // Найти число непроиндексированных узлов
        std::vector<int> indexedNodes;
        for (int i = 0; i < grid->outerNodes.size(); i++)
        {
            const auto node = grid->outerNodes[i];
            if (node->side == OuterNodeSide::Right)
            {
                continue;
            }
            if (node->side == OuterNodeSide::Top || node->side == OuterNodeSide::Bottom)
            {
                unindexedOuterNodes.emplace(grid->getNodeIndex(node));
                continue;
            }

            indexedNodes.push_back(grid->getNodeIndex(node));

            std::shared_ptr<Node> currentNode = node->parent.lock();
            while (1)
            {
                indexedNodes.push_back(grid->getNodeIndex(currentNode));

                auto currentNodeOuter = std::dynamic_pointer_cast<OuterNode>(currentNode);
                auto currentNodeInner = std::dynamic_pointer_cast<InnerNode>(currentNode);
                if (currentNodeOuter != nullptr)
                {
                    break;
                }
                currentNode = currentNodeInner->right.lock();
            }
        }
        // Учесть непроиндексированные внешние узлы, загрузить их в решение
        for (auto nodeIdx : unindexedOuterNodes)
        {
            T[iter][nodeIdx] = T[iter - 1][nodeIdx];
        }

        std::vector<int> indexedNodesRind(N, -1);
        for (int i = 0; i < indexedNodes.size(); i++)
        {
            indexedNodesRind[indexedNodes[i]] = i;
        }

        // Инициализировать трёхдиагональную матрицу
        TridiagonalMatrix matrix(indexedNodes.size());
        std::vector<double> vec(indexedNodes.size());

        // Составить СЛАУ для проиндексированных узлов
        for (int nodeIndex = 0; nodeIndex < grid->innerNodes.size(); nodeIndex++)
        {
            const auto node = grid->innerNodes[nodeIndex];

            const auto leftNode = node->left.lock();
            const auto rightNode = node->right.lock();
            const auto topNode = node->top.lock();
            const auto bottomNode = node->bottom.lock();

            const int leftNodeIndex = grid->getNodeIndex(leftNode);
            const int rightNodeIndex = grid->getNodeIndex(rightNode);
            const int topNodeIndex = grid->getNodeIndex(topNode);
            const int bottomNodeIndex = grid->getNodeIndex(bottomNode);

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

            double w_top = T[iter - 1][topNodeIndex];
            double w_bottom = T[iter - 1][bottomNodeIndex];

            int nodeIndexMat = indexedNodesRind[nodeIndex];
            int leftNodeIndexMat = indexedNodesRind[leftNodeIndex];
            int rightNodeIndexMat = indexedNodesRind[rightNodeIndex];

            matrix.set(nodeIndexMat, nodeIndexMat, -K * Kx * (mu_x + 1) - K * Ky * (mu_y + 1) - 1);
            matrix.set(nodeIndexMat, leftNodeIndexMat, K * Kx * mu_x);
            matrix.set(nodeIndexMat, rightNodeIndexMat, K * Kx);
            vec[nodeIndexMat] = -T[iter - 1][nodeIndex] - K * Ky * (mu_y * w_bottom + w_top);
        }
        for (int i = 0; i < grid->outerNodes.size(); i++)
        {
            const auto node = grid->outerNodes[i];
            if (node->side == OuterNodeSide::Bottom || node->side == OuterNodeSide::Top)
            {
                continue;
            }
            const auto parentNode = node->parent;
            const auto nodeIndex = grid->getNodeIndex(node);
            const auto parentNodeIndex = grid->getNodeIndex(node->parent.lock());
            const auto nodeIndexMat = indexedNodesRind[nodeIndex];
            const auto parentNodeIndexMat = indexedNodesRind[parentNodeIndex];

            const auto borderType = node->type;
            const auto relevantStep = ((node->side == OuterNodeSide::Left || node->side == OuterNodeSide::Right) ? grid->xStep : grid->yStep);
            assert(node->muValue > 0);
            assert(node->muValue <= 1);
            assert(nodeIndexMat != parentNodeIndexMat);
            if (borderType == BorderType::constFlow)
            {
                matrix.set(nodeIndexMat, parentNodeIndexMat, -1);
                matrix.set(nodeIndexMat, nodeIndexMat, 1);
                vec[nodeIndexMat] = node->borderValue * node->muValue * relevantStep;
            }
            else if (borderType == BorderType::constTemperature)
            {
                matrix.set(nodeIndexMat, nodeIndexMat, 1);
                vec[nodeIndexMat] = node->borderValue;
            }
            else if (borderType == BorderType::convection)
            {
                matrix.set(nodeIndexMat, nodeIndexMat, 1 - node->borderValue * node->muValue * relevantStep);
                matrix.set(nodeIndexMat, parentNodeIndexMat, -1);
                vec[nodeIndexMat] = 0;
            }
        }

        // Решить СЛАУ и загрузить результат в T
        auto solution = matrix.solve(vec);
        for (int i = 0; i < solution.size(); i++)
        {
            T[iter][indexedNodes[i]] = solution[i];
        }

        /** Вертикальный проход */
        unindexedOuterNodes.clear();
        // Найти число непроиндексированных узлов
        indexedNodes.clear();
        for (int i = 0; i < grid->outerNodes.size(); i++)
        {
            const auto node = grid->outerNodes[i];
            if (node->side == OuterNodeSide::Top)
            {
                continue;
            }
            if (node->side == OuterNodeSide::Left || node->side == OuterNodeSide::Right)
            {
                unindexedOuterNodes.emplace(grid->getNodeIndex(node));
                continue;
            }

            indexedNodes.push_back(grid->getNodeIndex(node));

            std::shared_ptr<Node> currentNode = node->parent.lock();
            while (1)
            {
                indexedNodes.push_back(grid->getNodeIndex(currentNode));

                auto currentNodeOuter = std::dynamic_pointer_cast<OuterNode>(currentNode);
                auto currentNodeInner = std::dynamic_pointer_cast<InnerNode>(currentNode);
                if (currentNodeOuter != nullptr)
                {
                    break;
                }
                currentNode = currentNodeInner->top.lock();
            }
        }

        indexedNodesRind = std::vector<int>(N, -1);
        for (int i = 0; i < indexedNodes.size(); i++)
        {
            indexedNodesRind[indexedNodes[i]] = i;
        }

        // Инициализировать трёхдиагональную матрицу
        matrix = TridiagonalMatrix(indexedNodes.size());
        vec = std::vector<double>(indexedNodes.size());

        // Составить СЛАУ для проиндексированных узлов
        for (int nodeIndex = 0; nodeIndex < grid->innerNodes.size(); nodeIndex++)
        {
            const auto node = grid->innerNodes[nodeIndex];

            const auto leftNode = node->left.lock();
            const auto rightNode = node->right.lock();
            const auto topNode = node->top.lock();
            const auto bottomNode = node->bottom.lock();

            const int leftNodeIndex = grid->getNodeIndex(leftNode);
            const int rightNodeIndex = grid->getNodeIndex(rightNode);
            const int topNodeIndex = grid->getNodeIndex(topNode);
            const int bottomNodeIndex = grid->getNodeIndex(bottomNode);

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

            double w_left = T[iter - 1][leftNodeIndex];
            double w_right = T[iter - 1][rightNodeIndex];

            int nodeIndexMat = indexedNodesRind[nodeIndex];
            int topNodeIndexMat = indexedNodesRind[topNodeIndex];
            int bottomNodeIndexMat = indexedNodesRind[bottomNodeIndex];

            matrix.set(nodeIndexMat, nodeIndexMat, -K * Kx * (mu_x + 1) - K * Ky * (mu_y + 1) - 1);
            matrix.set(nodeIndexMat, topNodeIndexMat, K * Kx * mu_x);
            matrix.set(nodeIndexMat, bottomNodeIndexMat, K * Kx);
            vec[nodeIndexMat] = -T[iter - 1][nodeIndex] - K * Ky * (mu_y * w_left + w_right);
        }
        for (int i = 0; i < grid->outerNodes.size(); i++)
        {
            const auto node = grid->outerNodes[i];
            if (node->side == OuterNodeSide::Left || node->side == OuterNodeSide::Right)
            {
                continue;
            }
            const auto parentNode = node->parent;
            const auto nodeIndex = grid->getNodeIndex(node);
            const auto parentNodeIndex = grid->getNodeIndex(node->parent.lock());
            const auto nodeIndexMat = indexedNodesRind[nodeIndex];
            const auto parentNodeIndexMat = indexedNodesRind[parentNodeIndex];

            const auto borderType = node->type;
            const auto relevantStep = ((node->side == OuterNodeSide::Left || node->side == OuterNodeSide::Right) ? grid->xStep : grid->yStep);
            assert(node->muValue > 0);
            assert(node->muValue <= 1);
            assert(nodeIndexMat != parentNodeIndexMat);
            if (borderType == BorderType::constFlow)
            {
                matrix.set(nodeIndexMat, parentNodeIndexMat, -1);
                matrix.set(nodeIndexMat, nodeIndexMat, 1);
                vec[nodeIndexMat] = node->borderValue * node->muValue * relevantStep;
            }
            else if (borderType == BorderType::constTemperature)
            {
                matrix.set(nodeIndexMat, nodeIndexMat, 1);
                vec[nodeIndexMat] = node->borderValue;
            }
            else if (borderType == BorderType::convection)
            {
                matrix.set(nodeIndexMat, nodeIndexMat, 1 - node->borderValue * node->muValue * relevantStep);
                matrix.set(nodeIndexMat, parentNodeIndexMat, -1);
                vec[nodeIndexMat] = 0;
            }
        }
        solution = matrix.solve(vec);
        for (int i = 0; i < solution.size(); i++)
        {
            T[iter][indexedNodes[i]] = solution[i];
        }
    }

    Solution sol;
    sol.T = T;
    sol.timeStep = dt;
    sol.grid = grid;
    return sol;
}
