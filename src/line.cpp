#include <algorithm>

#include "line.hpp"

std::vector<double> Line::xIntersections(double yLevel) const
{
    if (start.y == end.y)
    {
        return {};
    }

    const double x = -(start.x - end.x) / (start.y - end.y) * (end.y) + end.x;
    if (std::min(start.x, end.x) <= x &&
        std::max(start.x, end.x) >= x &&
        std::min(start.y, end.y) <= yLevel &&
        std::max(start.y, end.y) >= yLevel)
    {
        return {x};
    }

    return {};
}

std::vector<double> Line::yIntersections(double xLevel) const
{
    if (start.x == end.x)
    {
        return {};
    }

    const double y = -(start.y - end.y) / (start.x - end.x) * (end.x) + end.y;
    if (std::min(start.x, end.x) <= xLevel &&
        std::max(start.x, end.x) >= xLevel &&
        std::min(start.y, end.y) <= y &&
        std::max(start.y, end.y) >= y)
    {
        return {y};
    }

    return {};
}

std::pair<double, double> Line::xBoundaries() const
{
    return std::make_pair(std::min(start.x, end.x), std::max(start.x, end.x));
}

std::pair<double, double> Line::yBoundaries() const
{
    return std::make_pair(std::min(start.y, end.y), std::max(start.y, end.y));
}
