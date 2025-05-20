#include "solvers/solvers.hpp"

#include <cmath>
#include <set>
#include <iostream>

#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include "tridiagonal_matrix.hpp"
#include "solvers/solvers.hpp"
#include "solvers/implicit_solver.hpp"
#include "solvers/implicit_fast_solver.hpp"
#include "solvers/laplasian_calculator.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::solvers);

Solution solveImplicit(const Model &model, const std::shared_ptr<Grid> grid, double initT, double dt, int numIters)
{
    ImplicitSolver solver(model, grid, initT, dt, numIters);
    return solver.solve();
}

Solution solveImplicitFast(const Model &model, const std::shared_ptr<Grid> grid, double initT, double dt, int numIters)
{
    ImplicitFastSolver solver(model, grid, initT, dt, numIters);
    return solver.solve();
}

Solution calculateError(const Model &model, const std::shared_ptr<Grid> grid, const std::vector<double> &reference)
{
    ErrorCalculator calculator(model, grid, reference);
    return calculator.solve();
}

PLSAPR_END_NAMESPACE(); // plastinka_sapr::solvers
