#include <fstream>
#include <iostream>

#include "grid.hpp"

void Grid::writeToFile(const std::string &fileName) const
{
    std::ofstream f(fileName);

    f << "x y\n";

    for (auto node : innerNodes)
    {
        f << node->coords.x << " " << node->coords.y << "\n";
    }
    for (auto node : outerNodes)
    {
        f << node->coords.x << " " << node->coords.y << "\n";
    }

    f.close();
}

void Grid::writeToFile(const std::string &innerNodesFileName, const std::string &outerNodesFileName) const
{
    {
        std::ofstream f(innerNodesFileName);

        f << "x y\n";

        for (auto node : innerNodes)
        {
            f << node->coords.x << " " << node->coords.y << "\n";
        }

        f.close();
    }
    {
        std::ofstream f(outerNodesFileName);

        f << "x y\n";

        for (auto node : outerNodes)
        {
            f << node->coords.x << " " << node->coords.y << "\n";
        }

        f.close();
    }
}

int Grid::getNodeIndex(const std::shared_ptr<const Node> &node) const
{
    auto iter1 = innerIndex.find(node.get());
    auto iter2 = outerIndex.find(node.get());
    if (iter1 != innerIndex.end())
    {
        return iter1->second;
    }
    if (iter2 != outerIndex.end())
    {
        return iter2->second + static_cast<int>(innerNodes.size());
    }
    throw std::runtime_error("getNodeIndex failed");
}

void Grid::validateIntegrity(void) const
{
    for (int i = 0; i < innerNodes.size(); i++)
    {
        const auto &node = innerNodes[i];
        const auto &left = node->left;
        const auto &right = node->right;
        const auto &top = node->top;
        const auto &bottom = node->bottom;

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
}

void Grid::makeNodeIndexes(void)
{
    innerIndex = std::map<const Node *, int>{};
    outerIndex = std::map<const Node *, int>{};
    for (int i = 0; i < innerNodes.size(); i++)
    {
        outerIndex.emplace(innerNodes[i].get(), i);
    }
    for (int i = 0; i < outerNodes.size(); i++)
    {
        innerIndex.emplace(outerNodes[i].get(), i);
    }
}
