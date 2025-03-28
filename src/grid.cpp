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
