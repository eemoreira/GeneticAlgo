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
    SwapMutation(double rate) : mutationRate(rate) {}
    SwapMutation(const SwapMutation &) = default;
    SwapMutation(SwapMutation &&) = default;
    SwapMutation &operator=(const SwapMutation &) = default;
    SwapMutation &operator=(SwapMutation &&) = default;
};

struct PolynomialMutation {
    double mutationRate;
    double n_m; // distribution index
    template<typename Individual>
        void operator()(Individual& individual) {
            for (auto& gene : individual.genes) {
                if (std::uniform_real_distribution<>(0.0, 1.0)(rng) < mutationRate) {
                    double u = std::uniform_real_distribution<>(0.0, 1.0)(rng);
                    double delta = u < 0.5 ? std::pow(2 * u, 1.0 / (n_m + 1)) - 1 : 1 - std::pow(2 * (1 - u), 1.0 / (n_m + 1));
                    gene += delta;
                    gene = std::min(1.0, std::max(0.0, gene)); // assuming gene is in [0, 1]
                }
            }
        }
    PolynomialMutation(double rate, double eta) : mutationRate(rate), n_m(eta) {}
    PolynomialMutation(const PolynomialMutation &) = default;
    PolynomialMutation(PolynomialMutation &&) = default;
    PolynomialMutation &operator=(const PolynomialMutation &) = default;
    PolynomialMutation &operator=(PolynomialMutation &&) = default;
};
