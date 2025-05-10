#include <fstream>
#include <iostream>
#include <cassert>

#include "grid.hpp"

Point Node::getCoords(const std::shared_ptr<const Node> node)
{
    const auto innerNodePtr = dynamic_cast<const InnerNode *>(node.get());
    const auto outerNodePtr = dynamic_cast<const OuterNode *>(node.get());
    if (innerNodePtr != nullptr)
    {
        return innerNodePtr->coords;
    }
    if (outerNodePtr != nullptr)
    {
        return outerNodePtr->coords;
    }
    throw std::runtime_error("error while getting point coords");
}

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

int Grid::getNodeIndex(const std::shared_ptr<const Node> node) const
{
    auto iter1 = innerIndex.find(node.get());
    auto iter2 = outerIndex.find(node.get());
    if (iter1 != innerIndex.end())
    {
        const auto result = iter1->second;
        assert(result < innerNodes.size() + outerNodes.size());
        return result;
    }
    if (iter2 != outerIndex.end())
    {
        const auto result = iter2->second + static_cast<int>(innerNodes.size());
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

const std::shared_ptr<const Node> Grid::getNodeByIndex(int index) const
{
    if (index < innerNodes.size())
    {
        const auto result = std::dynamic_pointer_cast<const Node>(innerNodes[index]);
        assert(result != nullptr);
        return result;
    }
    if (index < getNodeCount())
    {
        const auto result = std::dynamic_pointer_cast<const Node>(outerNodes[index - static_cast<int>(innerNodes.size())]);
        assert(result != nullptr);
        return result;
    }
    throw std::runtime_error("error while getting node by index");
}

int Grid::getNodeCount(void) const
{
    return innerNodes.size() + outerNodes.size();
}
