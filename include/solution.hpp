#pragma once

#include <vector>
#include <memory>

#include "grid.hpp"

class Solution
{
public:
    std::vector<std::vector<double>> T; // Значения фазовой переменной
    std::shared_ptr<Grid> grid;
    double initialT;
    double timeStep;

    void printToFile(const std::string &fileName) const;
};
