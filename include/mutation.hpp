#pragma once

#include "rand.hpp"

struct BitFlipMutation {
    double mutationRate;
    template<typename Individual>
        void operator()(Individual& individual) {
            for (auto& gene : individual.genes) {
                if (std::uniform_real_distribution<>(0.0, 1.0)(rng) < mutationRate) {
                    gene = !gene;
                }
            }
        }
    BitFlipMutation(double rate) : mutationRate(rate) {}
    BitFlipMutation(const BitFlipMutation &) = default;
    BitFlipMutation(BitFlipMutation &&) = default;
    BitFlipMutation &operator=(const BitFlipMutation &) = default;
    BitFlipMutation &operator=(BitFlipMutation &&) = default;
    BitFlipMutation() = default;
};

struct SwapMutation {
    double mutationRate;
    template<typename Individual>
        void operator()(Individual& individual) {
            if (std::uniform_real_distribution<>(0.0, 1.0)(rng) < mutationRate) {
                size_t size = individual.genes.size();
                size_t i = rng() % size;
                size_t j = rng() % size;
                std::swap(individual.genes[i], individual.genes[j]);
            }
        }
    SwapMutation(const SwapMutation &) = default;
    SwapMutation(SwapMutation &&) = default;
    SwapMutation &operator=(const SwapMutation &) = default;
    SwapMutation &operator=(SwapMutation &&) = default;
};
