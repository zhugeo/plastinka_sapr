#include <iostream>
#include <memory>

#include "border.hpp"
#include "model.hpp"
#include "solvers.hpp"
#include "arc.hpp"
#include "line.hpp"
#include "grid_generator.hpp"

int main()
{
    Model model;
    model.borders.push_back(Border(
        std::make_shared<Line>(Point(0, 0), Point(0, 400)),
        BorderType::constTemperature, 100));
    model.borders.push_back(Border(
        std::make_shared<Line>(Point(0, 0), Point(500, 0)),
        BorderType::convection, 0.1));
    model.borders.push_back(Border(
        std::make_shared<Line>(Point(0, 400), Point(350, 400)),
        BorderType::constFlow, 0));
    model.borders.push_back(Border(
        std::make_shared<Line>(Point(500, 0), Point(500, 250)),
        BorderType::constTemperature, 100));

    model.borders.push_back(Border(
        std::make_shared<Arc>(Point(350, 250), 150, 0, 90),
        BorderType::constFlow, 0));
    model.borders.push_back(Border(
        std::make_shared<Arc>(Point(155, 255), 50, 0, 360),
        BorderType::constFlow, 0));

    auto g = generateGrid(model, 20, 20);
    g->writeToFile("innerNodes.csv", "outerNodes.csv");

    auto grid_pointer = std::shared_ptr<Grid>{std::move(g)};

    model.lambda = 100;
    model.c = 1;
    model.rho = 1;

    const auto solution = solveImplicit(model, grid_pointer, 1, 100, 0);
    solution.printToFile("solutionImplicit.csv");
}
