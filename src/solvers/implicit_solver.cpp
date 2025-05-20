#include "solvers/solvers.hpp"

#include <cmath>
#include <set>
#include <iostream>

#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include "tridiagonal_matrix.hpp"
#include "solvers/implicit_solver.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

std::vector<double> ImplicitSolver::solveStep(const std::vector<double> &prevT) const
{
    const int N = grid->getNodeCount();
    std::vector<double> thisT(N, 0);

    Eigen::SparseMatrix<double> matrix(N, N);
    Eigen::VectorXd vector = Eigen::VectorXd::Zero(N);

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

        matrix.insert(nodeIndex, nodeIndex) = (-K * Kx * (muX + 1) - K * Ky * (muY + 1) - 1);
        assert(!std::isnan(matrix.coeff(nodeIndex, nodeIndex)));
        matrix.insert(nodeIndex, leftNodeIndex) = K * Kx * muX;
        matrix.insert(nodeIndex, rightNodeIndex) = K * Kx;
        matrix.insert(nodeIndex, topNodeIndex) = K * Ky;
        matrix.insert(nodeIndex, bottomNodeIndex) = K * Ky * muY;
        vector(nodeIndex) = -prevT[nodeIndex];
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

    for (int i = 0; i < N; i++)
    {
        thisT[i] = solution(i);
    }
    return thisT;
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
