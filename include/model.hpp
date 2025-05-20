#pragma once

#include <vector>

#include "border.hpp"
#include "base.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr)

class Model
{
public:
    std::vector<Border> borders;
    double lambda, c, rho;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr
