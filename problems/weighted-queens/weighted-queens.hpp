#pragma once

#include <cstdint>
#include <vector>
#include <cmath>

struct WeightedQueensFitnessFunction {
    double penalty;
    std::vector<std::vector<double>> weights;
    double operator()(const std::vector<uint32_t>& genes) {
        uint32_t n = genes.size();
        uint32_t badPairs = 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                if (abs((int)genes[i] - (int)genes[j]) == abs((int)i - (int)j)) {
                    badPairs += 1;
                }
            }
        }

        double weightSum = 0;
        for (size_t i = 0; i < genes.size(); i++) {
            weightSum += weights[i][genes[i]];
        }
        return weightSum - penalty * badPairs;
    }
    WeightedQueensFitnessFunction(double p, const std::vector<std::vector<double>>& w) : penalty(p), weights(w) {}
};
