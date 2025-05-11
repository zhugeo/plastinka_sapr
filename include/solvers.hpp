#pragma once

#include <memory>

#include "grid.hpp"
#include "model.hpp"
#include "solution.hpp"

Solution solveImplicit(const Model &model, const std::shared_ptr<Grid> grid, double dt, int numIters, double initT);

Solution solveImplicitFast(const Model &model, const std::shared_ptr<Grid> grid, double dt, int numIters, double initT);
