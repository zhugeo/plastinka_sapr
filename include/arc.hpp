#pragma once

#include <algorithm>

#include "base.hpp"
#include "curve.hpp"

class Arc : public Curve
{
public:
    /** Угол равен нулю при направлении (0, 1); дальше по часовой стрелке */
    Arc(Point center, double R, double startAngle = 0.0, double stopAngle = 360.0) : center(center),
                                                                                     R(R),
                                                                                     startAngle(std::min(startAngle, stopAngle)),
                                                                                     stopAngle(std::max(startAngle, stopAngle)) {}
    virtual std::vector<double> xIntersections(double yLevel) const override;
    virtual std::vector<double> yIntersections(double xLevel) const override;

    virtual std::pair<double, double> xBoundaries() const override;
    virtual std::pair<double, double> yBoundaries() const override;
    virtual std::pair<double, double> getNormal(Point p) const override;

private:
    Point center;
    double R;
    double startAngle, stopAngle;
};
