#pragma once

#include <vector>
#include <memory>

class Curve
{
public:
    virtual std::vector<double> xIntersections(double yLevel) const = 0;
    virtual std::vector<double> yIntersections(double xLevel) const = 0;
    virtual std::pair<double, double> xBoundaries() const = 0;
    virtual std::pair<double, double> yBoundaries() const = 0;
    virtual ~Curve() {};
};

typedef std::shared_ptr<Curve> CurvePtr;
