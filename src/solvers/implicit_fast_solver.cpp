#include "solvers/implicit_fast_solver.hpp"

#include <cmath>
#include <set>
#include <iostream>
#include <cassert>

#include "matrices/tridiagonal_matrix.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

std::vector<double> ImplicitFastSolver::solveStep(const std::vector<double> &prevT) const
{
    const int N = grid->getNodeCount();
    std::vector<double> thisT(N, 0);

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
        thisT[nodeIdx] = prevT[nodeIdx];
    }

    std::vector<int> indexedNodesRind(N, -1);
    for (int i = 0; i < indexedNodes.size(); i++)
    {
        indexedNodesRind[indexedNodes[i]] = i;
    }

    // Инициализировать трёхдиагональную матрицу
    matrices::TridiagonalMatrix matrix(indexedNodes.size());
    std::vector<double> vec(indexedNodes.size());

    // Составить СЛАУ для проиндексированных узлов
    const auto &writeInnerNodeEquationHorizontal = [&](std::shared_ptr<InnerNode> node, double muX, double muY)
    {
        const auto nodeIndex = grid->getNodeIndex(node);

        const auto leftNode = node->left.lock();
        const auto rightNode = node->right.lock();
        const auto topNode = node->top.lock();
        const auto bottomNode = node->bottom.lock();

        const int leftNodeIndex = grid->getNodeIndex(leftNode);
        const int rightNodeIndex = grid->getNodeIndex(rightNode);
        const int topNodeIndex = grid->getNodeIndex(topNode);
        const int bottomNodeIndex = grid->getNodeIndex(bottomNode);

        double K = 2 * model.lambda * dt / model.c / model.rho;
        double Kx = 1 / (muX * (muX + 1) * grid->xStep * grid->xStep);
        double Ky = 1 / (muY * (muY + 1) * grid->yStep * grid->yStep);

        double w_top = prevT[topNodeIndex];
        double w_bottom = prevT[bottomNodeIndex];

        int nodeIndexMat = indexedNodesRind[nodeIndex];
        int leftNodeIndexMat = indexedNodesRind[leftNodeIndex];
        int rightNodeIndexMat = indexedNodesRind[rightNodeIndex];

        matrix.set(nodeIndexMat, nodeIndexMat, -K * Kx * (muX + 1) - K * Ky * (muY + 1) - 1);
        matrix.set(nodeIndexMat, leftNodeIndexMat, K * Kx * muX);
        matrix.set(nodeIndexMat, rightNodeIndexMat, K * Kx);
        vec[nodeIndexMat] = -prevT[nodeIndex] - K * Ky * (muY * w_bottom + w_top);
    };
    scanInnerNodes(writeInnerNodeEquationHorizontal);

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
        thisT[indexedNodes[i]] = solution[i];
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
    matrix = matrices::TridiagonalMatrix(indexedNodes.size());
    vec = std::vector<double>(indexedNodes.size());

    // Составить СЛАУ для проиндексированных узлов
    const auto &writeInnerNodeEquationVertical = [&](std::shared_ptr<InnerNode> node, double muX, double muY)
    {
        const auto nodeIndex = grid->getNodeIndex(node);

        const auto leftNode = node->left.lock();
        const auto rightNode = node->right.lock();
        const auto topNode = node->top.lock();
        const auto bottomNode = node->bottom.lock();

        const int leftNodeIndex = grid->getNodeIndex(leftNode);
        const int rightNodeIndex = grid->getNodeIndex(rightNode);
        const int topNodeIndex = grid->getNodeIndex(topNode);
        const int bottomNodeIndex = grid->getNodeIndex(bottomNode);

        double K = 2 * model.lambda * dt / model.c / model.rho;
        double Kx = 1 / (muX * (muX + 1) * grid->xStep * grid->xStep);
        double Ky = 1 / (muY * (muY + 1) * grid->yStep * grid->yStep);

        double w_left = prevT[leftNodeIndex];
        double w_right = prevT[rightNodeIndex];

        int nodeIndexMat = indexedNodesRind[nodeIndex];
        int topNodeIndexMat = indexedNodesRind[topNodeIndex];
        int bottomNodeIndexMat = indexedNodesRind[bottomNodeIndex];

        matrix.set(nodeIndexMat, nodeIndexMat, -K * Kx * (muX + 1) - K * Ky * (muY + 1) - 1);
        matrix.set(nodeIndexMat, topNodeIndexMat, K * Kx * muX);
        matrix.set(nodeIndexMat, bottomNodeIndexMat, K * Kx);
        vec[nodeIndexMat] = -prevT[nodeIndex] - K * Ky * (muY * w_left + w_right);
    };
    scanInnerNodes(writeInnerNodeEquationVertical);

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
        thisT[indexedNodes[i]] = solution[i];
    }
    return thisT;
}

Solution ImplicitFastSolver::solve()
{
    fillInitialStep();
    for (int iter = 1; iter < numIters; iter++)
    {
        T.push_back(solveStep(T.back()));
    }

    Solution s;
    s.grid = grid;
    s.initialT = initT;
    s.T = T;
    s.timeStep = dt;
    return s;
}

PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
