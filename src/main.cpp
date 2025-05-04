#include <iostream>
#include <memory>

#include "border.hpp"
#include "model.hpp"
#include "solvers.hpp"
#include "arc.hpp"
#include "line.hpp"
#include "generate_grid.hpp"

int main()
{
    Model model;
    model.borders.push_back(Border(
        std::make_shared<Arc>(Point(10, 10), 10, 0, 90),
        BorderType::constTemperature, 0));
    model.borders.push_back(Border(
        std::make_shared<Arc>(Point(10, 10), 10, 180, 270),
        BorderType::constTemperature, 0));
    model.borders.push_back(Border(
        std::make_shared<Line>(Point(0, 10), Point(0, 20)),
        BorderType::constTemperature, 0));
    model.borders.push_back(Border(
        std::make_shared<Line>(Point(10, 20), Point(0, 20)),
        BorderType::constTemperature, 0));
    model.borders.push_back(Border(
        std::make_shared<Line>(Point(10, 0), Point(20, 0)),
        BorderType::constTemperature, 0));
    model.borders.push_back(Border(
        std::make_shared<Line>(Point(20, 0), Point(20, 10)),
        BorderType::constTemperature, 0));

    Grid g = generateGrid(model, 2, 2);
    g.writeToFile("innerNodes.csv", "outerNodes.csv");

    solveImplicit(model, g, 0.1, 100, 100.0);
}
