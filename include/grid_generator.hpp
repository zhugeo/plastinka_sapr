#pragma once

#include <map>

#include "grid.hpp"
#include "model.hpp"

class GridGenerator
{
private:
    const Model &model;
    double xMax, xMin, yMax, yMin;
    double &xStep, &yStep;
    std::vector<double> xSlices, ySlices;

    std::map<std::pair<int, int>, std::shared_ptr<InnerNode>> innerNodes;

    std::map<std::pair<int, int>, std::shared_ptr<OuterNode>> outerHorizontalNodes;
    std::map<std::pair<int, int>, std::shared_ptr<OuterNode>> outerVerticalNodes;

    std::vector<std::shared_ptr<OuterNode>> outerNodes;

public:
    GridGenerator(const Model &model, double xStep, double yStep)
        : model(model),
          xStep(xStep),
          yStep(yStep)
    {
    }
    void generateGrid(void);

private:
    void calculateModelDimensions(void);
    void makeSlices(void);
    void bottomToUpScan(void);
    void leftToRightScan(void);
    void connectNodes(void);

    std::weak_ptr<Node> findNodeByCoords(int xSliceIndex,
                                         int ySliceIndex,
                                         bool includeHorizontalOuterNodes = false,
                                         bool includeVerticalOuterNodes = false) const;

    friend class GridExporter;
};

class GridExporter
{
private:
    const GridGenerator &generator;
    std::unique_ptr<Grid> grid;

public:
    GridExporter(const GridGenerator &generator) : generator(generator)
    {
    }
    std::unique_ptr<Grid> exportGrid(void);

private:
    void loadNodesToGrid(void);
    void validateGridIntegrity(void) const;
};

std::unique_ptr<Grid> generateGrid(const Model &model, double xStep, double yStep);
