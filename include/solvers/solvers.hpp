#pragma once

#include <memory>

#include "grid.hpp"
#include "model.hpp"
#include "solution.hpp"
#include "base.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

Solution solveExplicit(const Model &model, const std::shared_ptr<Grid> grid, double initT, double dt, int numIters);

Solution solveImplicit(const Model &model, const std::shared_ptr<Grid> grid, double initT, double dt, int numIters);

Solution solveImplicitFast(const Model &model, const std::shared_ptr<Grid> grid, double initT, double dt, int numIters);

Solution calculateLaplacian(const Model &model, const std::shared_ptr<Grid> grid, const std::vector<double> &reference);

PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
