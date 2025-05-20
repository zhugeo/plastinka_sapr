#include "solvers/solvers.hpp"

#include <cmath>
#include <set>
#include <iostream>

#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include "tridiagonal_matrix.hpp"
#include "solvers/laplacian_calculator.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

std::vector<double> LaplacianCalculator::solveStep(const std::vector<double> &prevT) const
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

        double K = 2 * model.lambda * 1 / model.c / model.rho;
        double Kx = 1 / (muX * (muX + 1) * grid->xStep * grid->xStep);
        double Ky = 1 / (muY * (muY + 1) * grid->yStep * grid->yStep);

        thisT[nodeIndex] = K * Kx * (muX * prevT[leftNodeIndex] - (muX + 1) * prevT[nodeIndex] + prevT[rightNodeIndex]) + K * Ky * (muY * prevT[bottomNodeIndex] - (muY + 1) * prevT[nodeIndex] + prevT[topNodeIndex]);
    };

    scanInnerNodes(writeInnerNodeEquation);

    return thisT;
}


Solution LaplacianCalculator::solve()
{
    std::vector<std::vector<double>> T;

    T.push_back(solveStep(reference));

    Solution s;
    s.grid = grid;
    s.T = T;
    return s;
}


PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
