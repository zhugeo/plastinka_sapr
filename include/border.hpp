#pragma once

#include "curve.hpp"

enum BorderType
{
    constTemperature,
    constFlow,
    convection,
};

class Border
{
public:
    Border(CurvePtr curve,
           BorderType type,
           double value) : curve(curve), type(type), value(value) {};
    CurvePtr curve;
    BorderType type;
    double value;
};
