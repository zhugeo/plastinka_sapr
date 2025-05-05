#include <vector>
#include <cmath>
#include <numbers>
#include <iostream>

#include "arc.hpp"

std::vector<double> Arc::xIntersections(double yLevel) const
{
    std::vector<double> candidatesX;
    if (yLevel < center.y - R || yLevel > center.y + R)
    {
        return {};
    }
    else if (yLevel == center.y - R || yLevel == center.y + R)
    {
        candidatesX.push_back(center.x);
    }
    else
    {
        candidatesX.push_back(center.x + std::sqrt(R * R - std::pow(yLevel - center.y, 2)));
        candidatesX.push_back(center.x - std::sqrt(R * R - std::pow(yLevel - center.y, 2)));
    }

    std::vector<double> answers;
    for (auto cand : candidatesX)
    {
        double angle = std::atan2(cand - center.x, yLevel - center.y);

        if (angle < 0)
        {
            angle += 2 * std::numbers::pi;
        }

        angle *= 180 / std::numbers::pi;

        if (angle >= startAngle - 0.1 && angle <= stopAngle + 0.1)
        {
            answers.push_back(cand);
        }
    }
    return answers;
}

std::vector<double> Arc::yIntersections(double xLevel) const
{
    std::vector<double> candidatesY;
    if (xLevel < center.x - R || xLevel > center.x + R)
    {
        return {};
    }
    else if (xLevel == center.x - R || xLevel == center.x + R)
    {
        candidatesY.push_back(center.y);
    }
    else
    {
        candidatesY.push_back(center.y + std::sqrt(R * R - std::pow(xLevel - center.x, 2)));
        candidatesY.push_back(center.y - std::sqrt(R * R - std::pow(xLevel - center.x, 2)));
    }

    std::vector<double> answers;
    for (auto cand : candidatesY)
    {
        double angle = std::atan2(xLevel - center.x, cand - center.y);

        if (angle < 0)
        {
            angle += 2 * std::numbers::pi;
        }

        angle *= 180 / std::numbers::pi;

        if (angle >= startAngle - 0.1 && angle <= stopAngle + 0.1)
        {
            answers.push_back(cand);
        }
    }
    return answers;
}

std::pair<double, double> Arc::xBoundaries() const
{
    return std::make_pair(center.x - R, center.x + R);
}

std::pair<double, double> Arc::yBoundaries() const
{
    return std::make_pair(center.y - R, center.y + R);
}

std::pair<double, double> Arc::getNormal(Point p) const
{
    double RLen = std::sqrt(std::pow(p.x - center.x, 2) + std::pow(p.x - center.x, 2));
    return std::make_pair((center.x - p.x) / RLen, (center.y - p.y) / RLen);
}
