#pragma once

#include <memory>
#include <functional>

#include "grid.hpp"
#include "model.hpp"
#include "solution.hpp"
#include "base.hpp"
#include "solvers/base_solver.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

class ImplicitFastSolver : public BaseSolver
{
private:
    double dt;
    int numIters;

    std::vector<double> solveStep(const std::vector<double> &prevT) const;

public:
    ImplicitFastSolver(const Model &model,
                       const std::shared_ptr<Grid> grid,
                       double initT,
                       double dt,
                       int numIters) : BaseSolver(model, grid, initT),
                                       dt(dt),
                                       numIters(numIters) {}
    virtual Solution solve() override;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
