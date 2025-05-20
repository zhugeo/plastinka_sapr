#pragma once

#include <memory>
#include <functional>

#include "grid.hpp"
#include "model.hpp"
#include "solution.hpp"
#include "base.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

class BaseSolver
{
protected:
    const Model &model;
    const std::shared_ptr<Grid> grid;
    double initT;
    std::vector<std::vector<double>> T;

    BaseSolver(const Model &model,
           const std::shared_ptr<Grid> grid,
           double initT) : model(model),
                           grid(grid),
                           initT(initT) {}

    /**
     * Обойти все внутренние узлы, для каждого из них вызвать fn.
     * Для каждого узла вычисляются muX и muY.
     * Длина верхнего смежного ребра в muX раз больше, чем нижнего
     * Длина правого смежного ребра в muY раз больше, чем левого
     */
    void scanInnerNodes(const std::function<void(std::shared_ptr<InnerNode> node, double muX, double muY)> &fn) const;

    void fillInitialStep(void);

public:
    virtual Solution solve() = 0;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
