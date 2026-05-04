#pragma once

#include <cmath>

struct CosFitnessFunction {
    double operator()(double x) {
        return cos(20 * x) - std::abs(x)/2 + (x*x*x)/4;
    }
};
