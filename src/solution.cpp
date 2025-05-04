#include "solution.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

void Solution::printToFile(const std::string &fileName) const
{
    std::stringstream out;

    out << "t x y T" << std::endl;
    for (const auto &timeLayer : this->T)
    {
    }
    std::ofstream f_out(fileName);
    f_out << out.str();
    f_out.close();
}
