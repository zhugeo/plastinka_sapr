#include "solvers/solvers.hpp"

#include <cmath>
#include <set>
#include <iostream>

#include "matrices/eigen_sparce_matrix.hpp"
#include "solvers/implicit_solver.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

std::vector<double> ImplicitSolver::solveStep(const std::vector<double> &prevT) const
{
    const auto N = grid->getNodeCount();

    auto matrix = matrices::SparceMatrixLU(N);
    auto vector = std::vector<double>(N, 0);

    const auto writeInnerNodeEquation = [&](std::shared_ptr<InnerNode> node, double muX, double muY)
    {
        const auto leftNode = node->left.lock();
        const auto rightNode = node->right.lock();
        const auto topNode = node->top.lock();
        const auto bottomNode = node->bottom.lock();

        const int nodeIndex = grid->getNodeIndex(node);
        const int leftNodeIndex = grid->getNodeIndex(node->left.lock());
        const int rightNodeIndex = grid->getNodeIndex(node->right.lock());
        const int topNodeIndex = grid->getNodeIndex(node->top.lock());
        const int bottomNodeIndex = grid->getNodeIndex(node->bottom.lock());

        double K = 2 * model.lambda * dt / model.c / model.rho;
        double Kx = 1 / (muX * (muX + 1) * grid->xStep * grid->xStep);
        double Ky = 1 / (muY * (muY + 1) * grid->yStep * grid->yStep);

        matrix.set(nodeIndex, nodeIndex, -K * Kx * (muX + 1) - K * Ky * (muY + 1) - 1);
        matrix.set(nodeIndex, leftNodeIndex, K * Kx * muX);
        matrix.set(nodeIndex, rightNodeIndex, K * Kx);
        matrix.set(nodeIndex, topNodeIndex, K * Ky);
        matrix.set(nodeIndex, bottomNodeIndex, K * Ky * muY);
        vector[nodeIndex] = -prevT[nodeIndex];
    };

    scanInnerNodes(writeInnerNodeEquation);

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
            matrix.set(nodeIndex, parentNodeIndex, -1);
            matrix.set(nodeIndex, nodeIndex, 1);
            vector[nodeIndex] = node->borderValue * node->muValue * relevantStep;
        }
        else if (borderType == BorderType::constTemperature)
        {
            matrix.set(nodeIndex, nodeIndex, 1);
            vector[nodeIndex] = node->borderValue;
        }
        else if (borderType == BorderType::convection)
        {
            matrix.set(nodeIndex, nodeIndex, 1 - node->borderValue * node->muValue * relevantStep);
            matrix.set(nodeIndex, parentNodeIndex, -1);
            vector[nodeIndex] = 0;
        }
    }

    return matrix.solve(vector);
}

Solution ImplicitSolver::solve()
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
