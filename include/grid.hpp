#pragma once

#include <vector>
#include <string>
#include <map>

#include "base.hpp"
#include "border.hpp"

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
    OuterNode(Point coords, OuterNodeSide side,
              BorderType type, double borderValue) : coords(coords),
                                                     side(side),
                                                     type(type),
                                                     borderValue(borderValue) {}
    Point coords;
    OuterNodeSide side;
    BorderType type;
    double borderValue;
};

class Grid
{
public:
    void writeToFile(const std::string &fileName) const;
    void writeToFile(const std::string &innerNodesFileName, const std::string &outerNodesFileName) const;

    int getNodeIndex(const std::shared_ptr<const Node> &) const;
    void makeNodeIndexes(void);
    void validateIntegrity(void) const;

    std::vector<std::shared_ptr<InnerNode>> innerNodes;
    std::vector<std::shared_ptr<OuterNode>> outerNodes;

    std::vector<double> xSlices;
    std::vector<double> ySlices;

    double xStep, yStep;

private:
    std::map<const Node *, int> innerIndex;
    std::map<const Node *, int> outerIndex;
};
