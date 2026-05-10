#pragma once

#include <cstdint>
#include <vector>

struct RadioFitnessFunction {
    double lambda;
    double penalty(double v0, double v1, double v2) const {
        return lambda * (v0 + v1 + v2);
    }
    std::pair<int32_t, int32_t> decode(const std::vector<uint8_t>& genes) const {
        uint32_t standard = 0;
        uint32_t luxury = 0;
        for (size_t i = 0; i < 5; i++) {
            standard = 2 * standard + genes[i];
        }
        for (size_t i = 5; i < 10; i++) {
            luxury = 2 * luxury + genes[i];
        }
        return {standard, luxury};
    }
    double operator()(const std::vector<uint8_t>& genes) {
        auto [standard, luxury] = decode(genes);
        int32_t v0 = std::max(0, standard - 24);
        int32_t v1 = std::max(0, luxury - 16);
        int32_t v2 = std::max(0, standard + 2 * luxury - 40);
        return std::max(0.0, standard * 30 + luxury * 40 - penalty(v0, v1, v2));
    }
    RadioFitnessFunction(const RadioFitnessFunction &) = default;
    RadioFitnessFunction(RadioFitnessFunction &&) = default;
    RadioFitnessFunction &operator=(const RadioFitnessFunction &) = default;
    RadioFitnessFunction &operator=(RadioFitnessFunction &&) = default;
    RadioFitnessFunction(double lambda) : lambda(lambda) {}
};
