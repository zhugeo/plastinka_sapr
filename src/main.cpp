#include <iostream>
#include <memory>

#include "border.hpp"
#include "model.hpp"
#include "solvers/solvers.hpp"
#include "curves/arc.hpp"
#include "curves/line.hpp"
#include "grid_generator.hpp"

using namespace plastinka_sapr;

int main()
{
    Model model;
    model.borders.push_back(Border(
        std::make_shared<curves::Line>(Point(0, 0), Point(0, 400)),
        BorderType::constTemperature, 100));
    model.borders.push_back(Border(
        std::make_shared<curves::Line>(Point(0, 0), Point(500, 0)),
        BorderType::convection, 1));
    model.borders.push_back(Border(
        std::make_shared<curves::Line>(Point(0, 400), Point(350, 400)),
        BorderType::convection, 1));
    model.borders.push_back(Border(
        std::make_shared<curves::Line>(Point(500, 0), Point(500, 250)),
        BorderType::constTemperature, 100));

    model.borders.push_back(Border(
        std::make_shared<curves::Arc>(Point(350, 250), 150, 0, 90),
        BorderType::convection, 1));
    model.borders.push_back(Border(
        std::make_shared<curves::Arc>(Point(155, 255), 50, 0, 360),
        BorderType::constFlow, 0));

    auto g = generateGrid(model, 10, 10);
    g->writeToFile("innerNodes.csv", "outerNodes.csv");

    auto grid_pointer = std::shared_ptr<Grid>{std::move(g)};

    model.lambda = 100;
    model.c = 1;
    model.rho = 1;

    const auto solution = solvers::solveExplicit(model, grid_pointer, 0.0, 1, 10);
    solution.printToFile("solutionImplicit.csv");

    // const auto sol2 = calculateError(model, grid_pointer, solution.T.back());
    // sol2.printToFile("errors.csv");
}
