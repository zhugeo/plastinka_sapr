#pragma once

#include <memory>
#include <functional>

#include "grid.hpp"
#include "model.hpp"
#include "solution.hpp"
#include "base.hpp"
#include "solvers/base_solver.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

class LaplacianCalculator : public BaseSolver
{
private:
    double dt;
    int numIters;
    const std::vector<double> &reference;

    std::vector<double> solveStep(const std::vector<double> &prevT) const;

public:
    LaplacianCalculator(const Model &model,
                    const std::shared_ptr<Grid> grid,
                    const std::vector<double> &reference) : BaseSolver(model, grid, 0),
                                                            reference(reference) {}
    virtual Solution solve() override;
};

Solution calculateError(const Model &model, const std::shared_ptr<Grid> grid, const std::vector<double> &reference);

PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
