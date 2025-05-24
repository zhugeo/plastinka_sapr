#pragma once

#include <memory>

#include "solvers/base_solver.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

class ImplicitSolver : public BaseSolver
{
private:
    double dt;
    int numIters;

    std::vector<double> solveStep(const std::vector<double> &prevT) const;

public:
    ImplicitSolver(const Model &model,
                   const std::shared_ptr<Grid> grid,
                   double initT,
                   double dt,
                   int numIters) : BaseSolver(model, grid, initT),
                                   dt(dt),
                                   numIters(numIters) {}
    virtual Solution solve() override;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
