#pragma once

#include <vector>
#include <memory>

#include "grid.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr);

class Solution
{
public:
    std::vector<std::vector<double>> T; // Значения фазовой переменной
    std::shared_ptr<Grid> grid;
    double initialT;
    double timeStep;

    void printToFile(const std::string &fileName) const;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr
