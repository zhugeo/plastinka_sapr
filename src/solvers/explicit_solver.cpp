#include "solvers/explicit_solver.hpp"

#include <cmath>
#include <set>
#include <cassert>

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

std::vector<double> ExplicitSolver::solveStep(const std::vector<double> &prevT) const
{
    const int N = grid->getNodeCount();
    std::vector<double> thisT(N, 0);

    const auto calculateNode = [&](std::shared_ptr<plastinka_sapr::InnerNode> node, double muX, double muY)
    {
        const double K = 2 * model.lambda;
        const double Kx = 1 / (muX * (muX + 1) * grid->xStep * grid->xStep);
        const double Ky = 1 / (muY * (muY + 1) * grid->yStep * grid->yStep);
        const int nodeIndex = grid->getNodeIndex(node);

        const double prevThisT = prevT[nodeIndex];
        const double prevTopT= prevT[grid->getNodeIndex(node->top.lock())];
        const double prevBottomT= prevT[grid->getNodeIndex(node->bottom.lock())];
        const double prevRightT= prevT[grid->getNodeIndex(node->right.lock())];
        const double prevLeftT= prevT[grid->getNodeIndex(node->left.lock())];

        thisT[nodeIndex] = prevThisT + K*Kx*muX*prevLeftT - (muX+1)*prevThisT + K*Kx*prevRightT + K*Ky*muY*prevBottomT
        - K*Ky*(muY+1)*prevThisT+K*Ky*prevTopT; };

    scanInnerNodes(calculateNode);

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
            thisT[nodeIndex] = node->borderValue * node->muValue * relevantStep + thisT[parentNodeIndex];
        }
        else if (borderType == BorderType::constTemperature)
        {
            thisT[nodeIndex] = node->borderValue;
        }
        else if (borderType == BorderType::convection)
        {
            thisT[nodeIndex] = thisT[parentNodeIndex] / (1 - node->borderValue * node->muValue * relevantStep);
        }
    }

    return thisT;
}

Solution ExplicitSolver::solve()
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
