#pragma once

#include <cstdint>
#include <vector>
#include <cmath>

struct NQueensFitnessFunction {
    uint32_t operator()(const std::vector<uint32_t>& genes) {
        uint32_t n = genes.size();
        uint32_t goodPairs = 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                if (abs((int)genes[i] - (int)genes[j]) != abs((int)i - (int)j)) {
                    goodPairs++;
                }
            }
        }
        return goodPairs;
    }
};
