#pragma once

#include <memory>
#include <functional>

#include "grid.hpp"
#include "model.hpp"
#include "solution.hpp"

class Solver
{
protected:
    const Model &model;
    const std::shared_ptr<Grid> grid;
    double initT;
    std::vector<std::vector<double>> T;

    Solver(const Model &model,
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

class ImplicitSolver : public Solver
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
                   int numIters) : Solver(model, grid, initT),
                                   dt(dt),
                                   numIters(numIters) {}
    virtual Solution solve() override;
};

class ImplicitFastSolver : public Solver
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
                       int numIters) : Solver(model, grid, initT),
                                       dt(dt),
                                       numIters(numIters) {}
    virtual Solution solve() override;
};

Solution solveImplicit(const Model &model, const std::shared_ptr<Grid> grid, double initT, double dt, int numIters);

Solution solveImplicitFast(const Model &model, const std::shared_ptr<Grid> grid, double initT, double dt, int numIters);
