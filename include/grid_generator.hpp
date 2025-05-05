#pragma once

#include <map>

#include "grid.hpp"
#include "model.hpp"

class GridGenerator
{
private:
    const Model &model;
    Grid grid;
    double xMax, xMin, yMax, yMin;
    double &xStep, &yStep;
    std::vector<double> &xSlices, &ySlices;

    std::map<std::pair<int, int>, std::shared_ptr<InnerNode>> innerNodes;

    std::map<std::pair<int, int>, std::shared_ptr<OuterNode>> outerHorizontalNodes;
    std::map<std::pair<int, int>, std::shared_ptr<OuterNode>> outerVerticalNodes;

public:
    GridGenerator(const Model &model, double xStep_, double yStep_)
        : model(model),
          xStep(grid.xStep),
          yStep(grid.yStep),
          xSlices(grid.xSlices),
          ySlices(grid.ySlices)
    {
        grid.xStep = xStep_;
        grid.yStep = yStep_;
    }
    Grid generateGrid(void);

private:
    void calculateModelDimensions(void);
    void makeSlices(void);
    void bottomToUpScan(void);
    void leftToRightScan(void);
    void connectNodes(void);
    void loadNodesToGrid(void);
    void validateGridIntegrity(void) const;

    std::weak_ptr<Node> findNodeByCoords(int xSliceIndex,
                                         int ySliceIndex,
                                         bool includeHorizontalOuterNodes = false,
                                         bool includeVerticalOuterNodes = false) const;
};

Grid generateGrid(const Model &model, double xStep, double yStep);
