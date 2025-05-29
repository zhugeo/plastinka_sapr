#pragma once

// https://stackoverflow.com/questions/3727862/is-there-any-way-to-make-visual-studio-stop-indenting-namespaces

#define PLSAPR_BEGIN_NAMESPACE(x) \
    namespace x                   \
    {

#define PLSAPR_END_NAMESPACE() }

#define PLSAPR_BEGIN_ANONYMOUS_NAMESPACE \
    namespace                            \
    {

class Point
{
public:
    Point(double x, double y) : x(x), y(y) {}

    double x;
    double y;
};
