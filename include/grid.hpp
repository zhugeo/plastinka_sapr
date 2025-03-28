#pragma once

#include <vector>
#include <string>

#include "base.hpp"

enum OuterNodeSide
{
    Left,
    Right,
    Top,
    Bottom
};

class Node
{
public:
    virtual ~Node() {}
};

class InnerNode : public Node
{
public:
    std::weak_ptr<Node> left, right, top, bottom;
    InnerNode(Point coords) : coords(coords) {}
    Point coords;
};

class OuterNode : public Node
{
public:
    std::weak_ptr<Node> left, right, top, bottom;
    OuterNode(Point coords, OuterNodeSide side) : coords(coords), side(side) {}
    Point coords;
    OuterNodeSide side;
};

class Grid
{
public:
    void writeToFile(const std::string &fileName) const;
    void writeToFile(const std::string &innerNodesFileName, const std::string &outerNodesFileName) const;

    std::vector<std::shared_ptr<InnerNode>> innerNodes;
    std::vector<std::shared_ptr<OuterNode>> outerNodes;

    std::vector<double> xSlices;
    std::vector<double> ySlices;
};
