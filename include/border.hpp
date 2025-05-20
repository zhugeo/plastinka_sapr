#pragma once

#include <memory>

#include "curves/curve.hpp"
#include "base.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr)

enum BorderType
{
    constTemperature,
    // Положительное направление потока - из тела
    constFlow,
    // Положительное направление потока - из тела
    convection,
};

class Border
{
public:
    Border(std::shared_ptr<const curves::Curve> curve,
           BorderType type,
           double value) : curve(curve), type(type), value(value) {};
    std::shared_ptr<const curves::Curve> curve;
    BorderType type;
    double value;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr
