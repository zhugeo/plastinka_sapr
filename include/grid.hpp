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
    static Point getCoords(const std::shared_ptr<const Node>);
    double muValue;
};

class InnerNode : public Node
{
public:
    std::weak_ptr<Node> left, right, top, bottom;
    InnerNode(Point coords) : coords(coords) { muValue = 1; }
    Point coords;
};

class OuterNode : public Node
{
public:
    std::weak_ptr<InnerNode> parent;
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

    int getNodeIndex(const std::shared_ptr<const Node>) const;
    const std::shared_ptr<const Node> getNodeByIndex(int index) const;

    int getNodeCount(void) const;

    std::vector<std::shared_ptr<InnerNode>> innerNodes;
    std::vector<std::shared_ptr<OuterNode>> outerNodes;

    std::vector<double> xSlices, ySlices;

    double xStep, yStep;

    friend class GridExporter;

private:
    void makeNodeIndexes(void);
    std::map<const Node *, int> innerIndex;
    std::map<const Node *, int> outerIndex;
};
