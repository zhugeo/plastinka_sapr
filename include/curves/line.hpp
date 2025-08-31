#pragma once

#include "curves/curve.hpp"
#include "base.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::curves);

class Line : public Curve
{
public:
    Line(Point start, Point end) : start(start), end(end) {}
    virtual std::vector<double> xIntersections(double yLevel) const override;
    virtual std::vector<double> yIntersections(double xLevel) const override;
    virtual std::pair<double, double> xBoundaries() const override;
    virtual std::pair<double, double> yBoundaries() const override;
    virtual std::pair<double, double> getNormal(Point p) const override;

private:
    Point start, end;
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr::curves
