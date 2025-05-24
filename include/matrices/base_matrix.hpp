#pragma once

#include <vector>

#include "base.hpp"

PLSAPR_BEGIN_NAMESPACE(plastinka_sapr::matrices);

class BaseMatrix
{
protected:
    int size;

public:
    virtual void set(int i, int j, double value) = 0;
    virtual int getSize() const final { return size; }
    virtual std::vector<double> solve(const std::vector<double> &x) const = 0;
    virtual ~BaseMatrix() {}

protected:
    BaseMatrix(int size) : size(size) {}
};

PLSAPR_END_NAMESPACE(); // plastinka_sapr::matrices
