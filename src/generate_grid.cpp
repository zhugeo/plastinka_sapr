#include <limits>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <map>

#include <iostream>

#include "generate_grid.hpp"

const double EPS = 1e-6;

bool approxEqual(double x, double y)
{
    return std::abs(x - y) < EPS;
}

Grid generateGrid(const Model &model, double xStep, double yStep)
{
    Grid answer;
    // Определить вернюю и нижнюю границы
    double xMax = std::numeric_limits<double>::min(),
           xMin = std::numeric_limits<double>::max(),
           yMax = std::numeric_limits<double>::min(),
           yMin = std::numeric_limits<double>::max();

    for (auto bord : model.borders)
    {
        auto xBound = bord.curve->xBoundaries();
        auto yBound = bord.curve->yBoundaries();
        xMax = std::max(xBound.second, xMax);
        yMax = std::max(yBound.second, yMax);
        xMin = std::min(xBound.first, xMin);
        yMin = std::min(yBound.first, yMin);
    }

    // Определим сечения по X и сечения по Y
    std::vector<double> xSlices, ySlices;
    for (double currentX = std::ceil(xMin / xStep) * xStep; currentX <= std::floor(xMax / xStep) * xStep; currentX += xStep)
    {
        xSlices.push_back(currentX);
    }
    for (double currentY = std::ceil(yMin / yStep) * yStep; currentY <= std::floor(yMax / yStep) * yStep; currentY += yStep)
    {
        ySlices.push_back(currentY);
    }

    answer.xSlices = xSlices;
    answer.ySlices = ySlices;

    // Создадим мапку для внутренних узлов
    std::map<std::pair<int, int>, std::shared_ptr<InnerNode>> innerNodes;

    // Создадим мапки для внешних узлов
    std::map<std::pair<int, int>, std::shared_ptr<OuterNode>> outerVerticalNodes;
    std::map<std::pair<int, int>, std::shared_ptr<OuterNode>> outerHorizontalNodes;

    // Проход "снизу вверх"
    for (int ySliceIndex = 0; ySliceIndex < ySlices.size(); ySliceIndex++)
    {
        auto currentY = ySlices[ySliceIndex];
        std::vector<double> intersections; // Массив X-координат точек пересечения
        for (auto bord : model.borders)
        {
            auto delta = bord.curve->xIntersections(currentY);
            for (auto in : delta)
            {
                intersections.push_back(in);
            }
        }
        std::sort(intersections.begin(), intersections.end());

        // Подчищаем "близко расположенные" точки
        for (int i = 0; i < intersections.size() - 1;)
        {
            if (approxEqual(intersections[i + 1], intersections[i]))
            {
                intersections.erase(intersections.begin() + i, intersections.begin() + i + 1);
            }
            else
            {
                i++;
            }
        }

        if (intersections.size() % 2 != 0)
        {
            throw std::runtime_error("Odd point count! Wrong boundary line");
        }

        // Создаём внутренние и внешние узлы

        int i = 0;

        for (int sliceIndex = 0; sliceIndex < xSlices.size() && i < intersections.size(); sliceIndex++)
        {
            auto currentX = xSlices[sliceIndex];

            // Проверить ред флаг
            if (i + 1 < intersections.size() && sliceIndex + 1 < xSlices.size())
            {
                if (intersections[i + 1] <= xSlices[sliceIndex + 1])
                {
                    throw std::runtime_error("Cant generate grid; try setting more fine grid step");
                }
            }

            if (approxEqual(intersections[i], currentX))
            {
                outerVerticalNodes.emplace(std::make_pair(sliceIndex, ySliceIndex),
                                           std::make_shared<OuterNode>(Point(intersections[i], currentY),
                                                                       (i % 2 == 0 ? OuterNodeSide::Bottom : OuterNodeSide::Top)));
                i++;
                continue;
            }
            else if (i % 2 != 0)
            {
                innerNodes.emplace(std::make_pair(sliceIndex, ySliceIndex),
                                   std::make_shared<InnerNode>(Point(currentX, currentY)));
            }
            if (intersections[i] < currentX + xStep)
            {
                outerVerticalNodes.emplace(std::make_pair(sliceIndex + (i % 2 == 0 ? 0 : +1), ySliceIndex),
                                           std::make_shared<OuterNode>(Point(intersections[i], currentY),
                                                                       (i % 2 == 0 ? OuterNodeSide::Bottom : OuterNodeSide::Top)));
                i++;
                continue;
            }
        }
    }

    // Проход "слева направо"
    for (int xSliceIndex = 0; xSliceIndex < xSlices.size(); xSliceIndex++)
    {
        auto currentX = xSlices[xSliceIndex];
        std::vector<double> intersections; // Массив Y-координат точек пересечения
        for (auto bord : model.borders)
        {
            auto delta = bord.curve->yIntersections(currentX);
            for (auto in : delta)
            {
                intersections.push_back(in);
                std::cout << in << std::endl;
            }
        }
        std::sort(intersections.begin(), intersections.end());
        std::cout << std::endl;

        // Подчищаем "близко расположенные" точки
        for (int i = 0; i < intersections.size() - 1;)
        {
            if (approxEqual(intersections[i + 1], intersections[i]))
            {
                intersections.erase(intersections.begin() + i, intersections.begin() + i + 1);
            }
            else
            {
                i++;
            }
        }

        if (intersections.size() % 2 != 0)
        {
            throw std::runtime_error("Odd point count! Wrong boundary line");
        }

        // Создаём внешние узлы и удаляем лишние внутренние

        int i = 0;

        for (int sliceIndex = 0; sliceIndex < ySlices.size() && i < intersections.size(); sliceIndex++)
        {
            auto currentY = ySlices[sliceIndex];

            // Проверить ред флаг
            if (i + 1 < intersections.size() && sliceIndex + 1 < ySlices.size())
            {
                if (intersections[i + 1] <= ySlices[sliceIndex + 1])
                {
                    throw std::runtime_error("Cant generate grid; try setting more fine grid step");
                }
            }

            if (approxEqual(intersections[i], currentY))
            {
                outerHorizontalNodes.emplace(std::make_pair(xSliceIndex, sliceIndex),
                                             std::make_shared<OuterNode>(Point(currentX, intersections[i]),
                                                                         (i % 2 == 0 ? OuterNodeSide::Left : OuterNodeSide::Right)));
                auto iter = innerNodes.find(std::make_pair(xSliceIndex, sliceIndex));
                if (iter != innerNodes.end())
                {
                    innerNodes.erase(iter);
                }
                i++;
                continue;
            }
            if (intersections[i] < currentY + yStep)
            {
                outerHorizontalNodes.emplace(std::make_pair(xSliceIndex, sliceIndex + (i % 2 == 0 ? 0 : +1)),
                                             std::make_shared<OuterNode>(Point(currentX, intersections[i]),
                                                                         (i % 2 == 0 ? OuterNodeSide::Left : OuterNodeSide::Right)));
                i++;
                continue;
            }
        }
    }

    // Обеспечим связность сетки
    for (auto node1 : innerNodes)
    {
        auto node = node1.second;
        answer.innerNodes.push_back(node);
        int xSliceIndex = node1.first.first;
        int ySliceIndex = node1.first.second;

        // Ищем нижнего соседа
        {

            auto iter1 = innerNodes.find(std::make_pair(xSliceIndex, ySliceIndex - 1));
            auto iter2 = outerHorizontalNodes.find(std::make_pair(xSliceIndex, ySliceIndex - 1));
            if (iter1 != innerNodes.end())
            {
                node->bottom = std::weak_ptr<Node>(iter1->second);
            }
            else if (iter2 != outerHorizontalNodes.end())
            {
                node->bottom = std::weak_ptr<Node>(iter2->second);
                iter2->second->top = std::weak_ptr<Node>(node);
                answer.outerNodes.push_back(iter2->second);
            }
            else
            {
                throw std::runtime_error("Error while connecting nodes");
            }
        }
        // Ищем верхнего соседа
        {

            auto iter1 = innerNodes.find(std::make_pair(xSliceIndex, ySliceIndex + 1));
            auto iter2 = outerHorizontalNodes.find(std::make_pair(xSliceIndex, ySliceIndex + 1));
            if (iter1 != innerNodes.end())
            {
                node->bottom = std::weak_ptr<Node>(iter1->second);
            }
            else if (iter2 != outerHorizontalNodes.end())
            {
                node->bottom = std::weak_ptr<Node>(iter2->second);
                iter2->second->bottom = std::weak_ptr<Node>(node);
                answer.outerNodes.push_back(iter2->second);
            }
            else
            {
                throw std::runtime_error("Error while connecting nodes");
            }
        }
        // Ищем правого соседа
        {

            auto iter1 = innerNodes.find(std::make_pair(xSliceIndex + 1, ySliceIndex));
            auto iter2 = outerVerticalNodes.find(std::make_pair(xSliceIndex + 1, ySliceIndex));
            if (iter1 != innerNodes.end())
            {
                node->bottom = std::weak_ptr<Node>(iter1->second);
            }
            else if (iter2 != outerVerticalNodes.end())
            {
                node->bottom = std::weak_ptr<Node>(iter2->second);
                iter2->second->left = std::weak_ptr<Node>(node);
                answer.outerNodes.push_back(iter2->second);
            }
            else
            {
                throw std::runtime_error("Error while connecting nodes");
            }
        }
        // Ищем левого соседа
        {

            auto iter1 = innerNodes.find(std::make_pair(xSliceIndex - 1, ySliceIndex));
            auto iter2 = outerVerticalNodes.find(std::make_pair(xSliceIndex - 1, ySliceIndex));
            if (iter1 != innerNodes.end())
            {
                node->bottom = std::weak_ptr<Node>(iter1->second);
            }
            else if (iter2 != outerVerticalNodes.end())
            {
                node->bottom = std::weak_ptr<Node>(iter2->second);
                iter2->second->right = std::weak_ptr<Node>(node);
                answer.outerNodes.push_back(iter2->second);
            }
            else
            {
                throw std::runtime_error("Error while connecting nodes");
            }
        }
    }

    return answer;
}
