#pragma once

#include <vector>
#include <memory>

#include "grid.hpp"

class Solution
{
public:
    std::vector<std::vector<double>> T; // Значения фазовой переменной
    Grid *grid;
    double initialT;

    void printToFile(const std::string &fileName) const;
};
