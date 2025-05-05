#include <fstream>
#include <iostream>
#include <cassert>

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
        const auto &result = iter1->second;
        assert(result < innerNodes.size() + outerNodes.size());
        return result;
    }
    if (iter2 != outerIndex.end())
    {
        const auto &result = iter2->second + static_cast<int>(innerNodes.size());
        assert(result < innerNodes.size() + outerNodes.size());
        return result;
    }
    throw std::runtime_error("getNodeIndex failed");
}

void Grid::makeNodeIndexes(void)
{
    innerIndex = std::map<const Node *, int>();
    outerIndex = std::map<const Node *, int>();
    for (int i = 0; i < innerNodes.size(); i++)
    {
        innerIndex.emplace(innerNodes[i].get(), i);
    }
    for (int i = 0; i < outerNodes.size(); i++)
    {
        outerIndex.emplace(outerNodes[i].get(), i);
    }
}
