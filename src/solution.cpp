#include "solution.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr);

void Solution::printToFile(const std::string &fileName) const
{
    std::stringstream out;

    out << "t x y T" << std::endl;
    for (int iterationNumber = 0; iterationNumber < T.size(); iterationNumber++)
    {
        for (int i = 0; i < grid->getNodeCount(); i++)
        {
            const auto node = grid->getNodeByIndex(i);
            const auto nodeCoords = Node::getCoords(node);
            out << iterationNumber * timeStep << " " << nodeCoords.x << " " << nodeCoords.y << " " << T[iterationNumber][i] << std::endl;
        }
    }
    std::ofstream f_out(fileName);
    f_out << out.str();
    f_out.close();
}

PLSAPR_END_NAMESPACE(); // plastinka_sapr
