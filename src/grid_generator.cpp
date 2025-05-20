#include "grid_generator.hpp"

#include <limits>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <cassert>

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr);

namespace
{
    const double EPS = 1e-6;

    bool approxEqual(double x, double y)
    {
        return std::abs(x - y) < EPS;
    }

    enum SliceDirection
    {
        parallelToX,
        parallelToY,
    };

    const std::vector<std::pair<double, Border>> sliceModel(const Model &model, double sliceLevel, SliceDirection direction)
    {
        std::vector<std::pair<double, Border>> intersections; // Массив X-координат точек пересечения
        for (auto border : model.borders)
        {
            std::vector<double> delta;
            if (direction == SliceDirection::parallelToX)
            {
                delta = border.curve->xIntersections(sliceLevel);
            }
            else
            {
                delta = border.curve->yIntersections(sliceLevel);
            }

            for (auto in : delta)
            {
                intersections.push_back(std::make_pair(in, border));
            }
        }

        std::sort(intersections.begin(), intersections.end(), [](auto a, auto b)
                  { return a.first < b.first; });

        // Подчищаем "близко расположенные" точки
        for (int i = 0; i < intersections.size() - 1;)
        {
            if (approxEqual(intersections[i + 1].first, intersections[i].first))
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

        return intersections;
    }
} // namespace

void GridGenerator::calculateModelDimensions(void)
{
    xMax = std::numeric_limits<double>::min();
    xMin = std::numeric_limits<double>::max();
    yMax = std::numeric_limits<double>::min();
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
}

void GridGenerator::makeSlices(void)
{
    for (double currentX = std::ceil(xMin / xStep) * xStep; currentX <= std::floor(xMax / xStep) * xStep; currentX += xStep)
    {
        xSlices.push_back(currentX);
    }
    for (double currentY = std::ceil(yMin / yStep) * yStep; currentY <= std::floor(yMax / yStep) * yStep; currentY += yStep)
    {
        ySlices.push_back(currentY);
    }
}

void GridGenerator::bottomToUpScan(void)
{
    // Проход "снизу вверх"
    for (int ySliceIndex = 0; ySliceIndex < ySlices.size(); ySliceIndex++)
    {
        double currentY = ySlices[ySliceIndex];
        const auto intersections = sliceModel(model, currentY, SliceDirection::parallelToX);

        // Создаём внутренние и внешние узлы

        int i = 0;

        for (int xSliceIndex = 0; xSliceIndex < xSlices.size() && i < intersections.size(); xSliceIndex++)
        {
            const auto &currentX = xSlices[xSliceIndex];

            // Проверить, является ли модель "плохой" для создания сетки
            if (i + 1 < intersections.size() && xSliceIndex + 1 < xSlices.size())
            {
                if (intersections[i + 1].first <= xSlices[xSliceIndex + 1])
                {
                    throw std::runtime_error("Cant generate grid; try setting more fine grid step");
                }
            }

            if (approxEqual(intersections[i].first, currentX))
            {
                Point point(intersections[i].first, currentY);
                auto border = intersections[i].second;
                auto borderValue = border.value;
                if (border.type != BorderType::constTemperature)
                {
                    auto normal = border.curve->getNormal(point);
                    if ((i % 2 == 0 && normal.second > 0) || (i % 2 != 0 && normal.second < 0))
                    {
                        normal.first = -normal.first;
                        normal.second = -normal.second;
                    }
                    borderValue *= normal.first;
                }
                outerHorizontalNodes.emplace(std::make_pair(xSliceIndex, ySliceIndex),
                                             std::make_shared<OuterNode>(point,
                                                                         (i % 2 == 0 ? OuterNodeSide::Left : OuterNodeSide::Right),
                                                                         border.type,
                                                                         borderValue));
                i++;
                continue;
            }
            else if (i % 2 != 0)
            {
                innerNodes.emplace(std::make_pair(xSliceIndex, ySliceIndex),
                                   std::make_shared<InnerNode>(Point(currentX, currentY)));
            }
            if (intersections[i].first < currentX + xStep)
            {
                Point point(intersections[i].first, currentY);
                auto border = intersections[i].second;
                auto borderValue = border.value;
                if (border.type != BorderType::constTemperature)
                {
                    auto normal = border.curve->getNormal(point);
                    if ((i % 2 == 0 && normal.second > 0) || (i % 2 != 0 && normal.second < 0))
                    {
                        normal.first = -normal.first;
                        normal.second = -normal.second;
                    }
                    borderValue *= normal.first;
                }
                outerHorizontalNodes.emplace(std::make_pair(xSliceIndex + (i % 2 == 0 ? 0 : +1), ySliceIndex),
                                             std::make_shared<OuterNode>(point,
                                                                         (i % 2 == 0 ? OuterNodeSide::Left : OuterNodeSide::Right),
                                                                         border.type,
                                                                         borderValue));
                i++;
                continue;
            }
        }
    }
}

void GridGenerator::leftToRightScan(void)
{
    // Проход "слева направо"
    for (int xSliceIndex = 0; xSliceIndex < xSlices.size(); xSliceIndex++)
    {
        auto currentX = xSlices[xSliceIndex];
        const auto intersections = sliceModel(model, currentX, SliceDirection::parallelToY);

        // Создаём внешние узлы и удаляем лишние внутренние

        int i = 0;

        for (int ySliceIndex = 0; ySliceIndex < ySlices.size() && i < intersections.size(); ySliceIndex++)
        {
            auto currentY = ySlices[ySliceIndex];

            // Проверить ред флаг
            if (i + 1 < intersections.size() && ySliceIndex + 1 < ySlices.size())
            {
                if (intersections[i + 1].first <= ySlices[ySliceIndex + 1])
                {
                    throw std::runtime_error("Cant generate grid; try setting more fine grid step");
                }
            }

            if (approxEqual(intersections[i].first, currentY))
            {
                Point point(currentX, intersections[i].first);
                auto border = intersections[i].second;
                auto borderValue = border.value;
                if (border.type != BorderType::constTemperature)
                {
                    auto normal = border.curve->getNormal(point);
                    if ((i % 2 == 0 && normal.second > 0) || (i % 2 != 0 && normal.second < 0))
                    {
                        normal.first = -normal.first;
                        normal.second = -normal.second;
                    }
                    borderValue *= normal.second;
                }
                const auto sliceCoords = std::make_pair(xSliceIndex, ySliceIndex);
                outerVerticalNodes.emplace(sliceCoords,
                                           std::make_shared<OuterNode>(point,
                                                                       (i % 2 == 0 ? OuterNodeSide::Bottom : OuterNodeSide::Top),
                                                                       border.type,
                                                                       borderValue));
                auto iter = innerNodes.find(sliceCoords);
                if (iter != innerNodes.end())
                {
                    innerNodes.erase(iter);
                }
                i++;
                continue;
            }
            if (intersections[i].first < currentY + yStep)
            {
                Point point(currentX, intersections[i].first);
                auto border = intersections[i].second;
                auto borderValue = border.value;
                if (border.type != BorderType::constTemperature)
                {
                    auto normal = border.curve->getNormal(point);
                    if ((i % 2 == 0 && normal.second > 0) || (i % 2 != 0 && normal.second < 0))
                    {
                        normal.first = -normal.first;
                        normal.second = -normal.second;
                    }
                    borderValue *= normal.second;
                }
                outerVerticalNodes.emplace(std::make_pair(xSliceIndex, ySliceIndex + (i % 2 == 0 ? 0 : +1)),
                                           std::make_shared<OuterNode>(point,
                                                                       (i % 2 == 0 ? OuterNodeSide::Bottom : OuterNodeSide::Top),
                                                                       border.type,
                                                                       borderValue));
                i++;
                continue;
            }
        }
    }
}

void GridGenerator::connectNodes(void)
{
    const auto connectNodeIfOuter = [&](std::shared_ptr<InnerNode> parentNode, std::weak_ptr<Node> node)
    {
        auto lockedNode = node.lock();
        const auto outerNodeCast = std::dynamic_pointer_cast<OuterNode>(lockedNode);
        if (outerNodeCast != nullptr)
        {
            outerNodes.push_back(outerNodeCast);
            outerNodeCast->parent = std::weak_ptr(parentNode);
            const auto nodeCoords = outerNodeCast->coords;
            const auto parentNodeCoords = parentNode->coords;
            if (approxEqual(nodeCoords.x, parentNodeCoords.x))
            {
                outerNodeCast->muValue = std::abs(nodeCoords.y - parentNodeCoords.y) / yStep;
            }
            else
            {
                outerNodeCast->muValue = std::abs(nodeCoords.x - parentNodeCoords.x) / xStep;
            }
        }
    };
    // Обеспечим связность сетки
    for (auto node1 : innerNodes)
    {
        auto node = node1.second;
        const auto xSliceIndex = node1.first.first;
        const auto ySliceIndex = node1.first.second;

        auto bottom = findNodeByCoords(xSliceIndex, ySliceIndex - 1, false, true);
        auto top = findNodeByCoords(xSliceIndex, ySliceIndex + 1, false, true);
        auto right = findNodeByCoords(xSliceIndex + 1, ySliceIndex, true, false);
        auto left = findNodeByCoords(xSliceIndex - 1, ySliceIndex, true, false);

        node->bottom = bottom;
        node->top = top;
        node->right = right;
        node->left = left;

        connectNodeIfOuter(node, bottom);
        connectNodeIfOuter(node, top);
        connectNodeIfOuter(node, right);
        connectNodeIfOuter(node, left);
    }
}

std::weak_ptr<Node> GridGenerator::findNodeByCoords(
    int xSliceIndex, int ySliceIndex,
    bool includeHorizontalOuterNodes, bool includeVerticalOuterNodes) const
{
    const auto searchTemplate = std::make_pair(xSliceIndex, ySliceIndex);
    const auto innerNodeIterator = innerNodes.find(searchTemplate);
    if (innerNodeIterator != innerNodes.end())
    {
        return std::weak_ptr<Node>(innerNodeIterator->second);
    }

    if (includeVerticalOuterNodes)
    {
        const auto outerVerticalNodeIterator = outerVerticalNodes.find(searchTemplate);
        if (outerVerticalNodeIterator != outerVerticalNodes.end())
        {
            return std::weak_ptr<Node>(outerVerticalNodeIterator->second);
        }
    }

    if (includeHorizontalOuterNodes)
    {
        const auto outerHorizontalNodeIterator = outerHorizontalNodes.find(searchTemplate);
        if (outerHorizontalNodeIterator != outerHorizontalNodes.end())
        {
            return std::weak_ptr<Node>(outerHorizontalNodeIterator->second);
        }
    }

    throw std::runtime_error("Node with such coords not found");
}

void GridExporter::loadNodesToGrid(void)
{
    for (auto i : generator.innerNodes)
    {
        const auto innerNode = i.second;
        grid->innerNodes.push_back(innerNode);
    }
    for (auto outerNode : generator.outerNodes)
    {
        grid->outerNodes.push_back(outerNode);
    }
}

void GridExporter::validateGridIntegrity(void) const
{
    for (int i = 0; i < grid->innerNodes.size(); i++)
    {
        const auto node = grid->innerNodes[i];

        const auto left = node->left;
        const auto right = node->right;
        const auto top = node->top;
        const auto bottom = node->bottom;

        if (left.expired())
        {
            throw std::runtime_error("Left neighbour is missing");
        }
        if (right.expired())
        {
            throw std::runtime_error("Right neighbour is missing");
        }
        if (top.expired())
        {
            throw std::runtime_error("Top neighbour is missing");
        }
        if (bottom.expired())
        {
            throw std::runtime_error("Bottom neighbour is missing");
        }
    }
    for (int i = 0; i < grid->outerNodes.size(); i++)
    {
        if (grid->outerNodes[i]->parent.expired())
        {
            throw std::runtime_error("Parent is missing");
        }
    }
    for (int i = 0; i < grid->outerNodes.size(); i++)
    {
        const auto node = grid->outerNodes[i];
        const auto parentNode = node->parent.lock();
        if (node->side == OuterNodeSide::Left)
        {
            assert(parentNode->left.lock() == node);
        }
        if (node->side == OuterNodeSide::Right)
        {
            assert(parentNode->right.lock() == node);
        }
        if (node->side == OuterNodeSide::Top)
        {
            assert(parentNode->top.lock() == node);
        }
        if (node->side == OuterNodeSide::Bottom)
        {
            assert(parentNode->bottom.lock() == node);
        }
    }
}

void GridGenerator::generateGrid(void)
{
    calculateModelDimensions();
    makeSlices();
    bottomToUpScan();
    leftToRightScan();
    connectNodes();
}

std::unique_ptr<Grid> GridExporter::exportGrid(void)
{
    grid = std::make_unique<Grid>();
    grid->xSlices = generator.xSlices;
    grid->ySlices = generator.ySlices;
    grid->xStep = generator.xStep;
    grid->yStep = generator.yStep;

    loadNodesToGrid();
    validateGridIntegrity();
    grid->makeNodeIndexes();
    return std::move(grid);
}

std::unique_ptr<Grid> generateGrid(const Model &model, double xStep, double yStep)
{
    auto generator = GridGenerator(model, xStep, yStep);
    generator.generateGrid();
    auto exporter = GridExporter(generator);
    return exporter.exportGrid();
}

PLSAPR_END_NAMESPACE(); // plastinka_sapr
