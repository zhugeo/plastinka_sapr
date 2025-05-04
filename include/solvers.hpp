#pragma once

#include "grid.hpp"
#include "model.hpp"
#include "solution.hpp"

Solution solveImplicit(const Model &model, const Grid &grid, double dt, int numIters, double initT);
