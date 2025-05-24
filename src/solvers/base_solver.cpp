#include "solvers/base_solver.hpp"

#include <cmath>
#include <set>

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

void BaseSolver::scanInnerNodes(const std::function<void(std::shared_ptr<InnerNode> node, double muX, double muY)> &fn) const
{
    for (int i = 0; i < grid->innerNodes.size(); i++)
    {
        auto node = grid->innerNodes[i];

        auto leftNode = node->left.lock();
        auto rightNode = node->right.lock();
        auto topNode = node->top.lock();
        auto bottomNode = node->bottom.lock();

        double muX = rightNode->muValue;
        if (muX == 1)
        {
            muX = 1 / leftNode->muValue;
        }

        double muY = topNode->muValue;
        if (muY == 1)
        {
            muY = 1 / bottomNode->muValue;
        }

        fn(node, muX, muY);
    }
}

void BaseSolver::fillInitialStep(void)
{
    T = std::vector<std::vector<double>>();
    int N = grid->getNodeCount();

    T.push_back(std::vector<double>(N, initT));
}

PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
