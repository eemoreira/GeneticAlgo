#pragma once

#include <vector>
#include <algorithm>
#include "rand.hpp"

struct UniformCrossover {
    double crossoverRate;
    template<typename Individual>
        std::pair<Individual, Individual> operator()(
            const Individual& parent1,
            const Individual& parent2
        ) {
            if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > crossoverRate) {
                return {parent1, parent2};
            }
            size_t size = parent1.genes.size();
            Individual child1, child2;
            child1.genes.resize(size);
            child2.genes.resize(size);

            for (size_t i = 0; i < size; ++i) {
                if (std::uniform_real_distribution(0.0, 1.0)(rng) < 0.5) {
                    child1.genes[i] = parent1.genes[i];
                    child2.genes[i] = parent2.genes[i];
                } else {
                    child1.genes[i] = parent2.genes[i];
                    child2.genes[i] = parent1.genes[i];
                }
            }

            return {child1, child2};
        }
    UniformCrossover(double rate) : crossoverRate(rate) {}
    UniformCrossover(const UniformCrossover &) = default;
    UniformCrossover(UniformCrossover &&) = default;
    UniformCrossover &operator=(const UniformCrossover &) = default;
    UniformCrossover &operator=(UniformCrossover &&) = default;
    UniformCrossover() = default;
};

struct OnePointCrossover {
    double crossoverRate;
    template<typename Individual>
        std::pair<Individual, Individual> operator()(
            const Individual& parent1,
            const Individual& parent2
        ) {
            if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > crossoverRate) {
                return {parent1, parent2};
            }
            size_t size = parent1.genes.size();
            size_t crossoverPoint = rng() % size;

            Individual child1, child2;
            child1.genes.resize(size);
            child2.genes.resize(size);

            for (size_t i = 0; i < crossoverPoint; ++i) {
                child1.genes[i] = parent1.genes[i];
                child2.genes[i] = parent2.genes[i];
            }
            for (size_t i = crossoverPoint; i < size; ++i) {
                child1.genes[i] = parent2.genes[i];
                child2.genes[i] = parent1.genes[i];
            }

            return {child1, child2};
        }
    OnePointCrossover(double rate) : crossoverRate(rate) {}
    OnePointCrossover(const OnePointCrossover &) = default;
    OnePointCrossover(OnePointCrossover &&) = default;
    OnePointCrossover &operator=(const OnePointCrossover &) = default;
    OnePointCrossover &operator=(OnePointCrossover &&) = default;
    OnePointCrossover() = default;
};

struct TwoPointCrossover {
    double crossoverRate;
    template<typename Individual>
        std::pair<Individual, Individual> operator()(
            const Individual& parent1,
            const Individual& parent2
        ) {
            if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > crossoverRate) {
                return {parent1, parent2};
            }
            size_t size = parent1.genes.size();
            size_t point1 = rng() % size;
            size_t point2 = rng() % size;

            if (point1 > point2) std::swap(point1, point2);

            Individual child1, child2;
            child1.genes.resize(size);
            child2.genes.resize(size);

            for (size_t i = 0; i < point1; ++i) {
                child1.genes[i] = parent1.genes[i];
                child2.genes[i] = parent2.genes[i];
            }
            for (size_t i = point1; i < point2; ++i) {
                child1.genes[i] = parent2.genes[i];
                child2.genes[i] = parent1.genes[i];
            }
            for (size_t i = point2; i < size; ++i) {
                child1.genes[i] = parent1.genes[i];
                child2.genes[i] = parent2.genes[i];
            }

            return {child1, child2};
        }
    TwoPointCrossover(double rate) : crossoverRate(rate) {}
    TwoPointCrossover(const TwoPointCrossover &) = default;
    TwoPointCrossover(TwoPointCrossover &&) = default;
    TwoPointCrossover &operator=(const TwoPointCrossover &) = default;
    TwoPointCrossover &operator=(TwoPointCrossover &&) = default;
    TwoPointCrossover() = default;
};

struct PMXCrossover {
    double crossoverRate;
    template<typename Individual>
        std::pair<Individual, Individual> operator()(
            const Individual& parent1,
            const Individual& parent2
        ) {
            if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > crossoverRate) {
                return {parent1, parent2};
            }

            size_t size = parent1.genes.size();
            size_t point1 = rng() % size;
            size_t point2 = rng() % size;

            if (point1 > point2) {
                std::swap(point1, point2);
            }

            Individual child1, child2;

            child1.genes.resize(size, -1);
            child2.genes.resize(size, -1);

            // Copy crossover segment
            for (size_t i = point1; i <= point2; ++i) {
                child1.genes[i] = parent2.genes[i];
                child2.genes[i] = parent1.genes[i];
            }

            auto pmx_fill = [&](auto& child,
                    const auto& parentA,
                    const auto& parentB) {

                for (size_t i = 0; i < size; ++i) {

                    if (i >= point1 && i <= point2)
                        continue;

                    auto gene = parentA.genes[i];

                    while (std::find(
                            child.genes.begin() + point1,
                            child.genes.begin() + point2 + 1,
                            gene
                    ) != child.genes.begin() + point2 + 1) {

                        auto it = std::find(
                                parentB.genes.begin() + point1,
                                parentB.genes.begin() + point2 + 1,
                                gene
                                );

                        size_t idx = std::distance(parentB.genes.begin(), it);

                        gene = parentA.genes[idx];
                    }

                    child.genes[i] = gene;
                }
            };

            pmx_fill(child1, parent1, parent2);
            pmx_fill(child2, parent2, parent1);

            return {child1, child2};
        }


    PMXCrossover(double rate) : crossoverRate(rate) {}
    PMXCrossover(const PMXCrossover &) = default;
    PMXCrossover(PMXCrossover &&) = default;
    PMXCrossover &operator=(const PMXCrossover &) = default;
    PMXCrossover &operator=(PMXCrossover &&) = default;
    PMXCrossover() = default;
};

struct CXCrossover {
    double crossoverRate;
    template<typename Individual>
        std::pair<Individual, Individual> operator()(
            const Individual& parent1,
            const Individual& parent2
        ) {
            if (std::uniform_real_distribution(0.0, 1.0)(rng) > crossoverRate) {
                return {parent1, parent2};
            }
            size_t size = parent1.genes.size();
            Individual child1, child2;
            child1.genes.resize(size);
            child2.genes.resize(size);

            std::vector<uint8_t> visited(size, false);
            for (size_t i = 0; i < size; ++i) {
                if (!visited[i]) {
                    size_t j = i;
                    do {
                        child1.genes[j] = parent1.genes[j];
                        child2.genes[j] = parent2.genes[j];
                        visited[j] = true;
                        j = std::find(parent1.genes.begin(), parent1.genes.end(), parent2.genes[j]) - parent1.genes.begin();
                    } while (j != i);
                }
            }

            return {child1, child2};
        }
    CXCrossover(double rate) : crossoverRate(rate) {}
    CXCrossover(const CXCrossover &) = default;
    CXCrossover(CXCrossover &&) = default;
    CXCrossover &operator=(const CXCrossover &) = default;
    CXCrossover &operator=(CXCrossover &&) = default;
    CXCrossover() = default;
};

struct SBXCrossover {
    double crossoverRate;
    double n_c; // distribution index
    template<typename Individual>
        std::pair<Individual, Individual> operator()(
            const Individual& parent1,
            const Individual& parent2
        ) {
            if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) > crossoverRate) {
                return {parent1, parent2};
            }
            size_t genes_size = parent1.genes.size();
            size_t fitness_size = parent1.fitness.size();
            Individual child1, child2;
            child1.genes.resize(genes_size);
            child2.genes.resize(genes_size);
            child1.fitness.resize(fitness_size);
            child2.fitness.resize(fitness_size);

            for (size_t i = 0; i < genes_size; ++i) {
                double u = std::uniform_real_distribution<>(0.0, 1.0)(rng);
                double beta = u < 0.5 ? std::pow(2 * u, 1.0 / (n_c + 1)) : std::pow(2 * (1 - u), -1.0 / (n_c + 1));
                child1.genes[i] = 0.5 * ((1 + beta) * parent1.genes[i] + (1 - beta) * parent2.genes[i]);
                child2.genes[i] = 0.5 * ((1 - beta) * parent1.genes[i] + (1 + beta) * parent2.genes[i]);
                child1.genes[i] = std::min(1.0, std::max(0.0, child1.genes[i]));
                child2.genes[i] = std::min(1.0, std::max(0.0, child2.genes[i]));
            }

            return {child1, child2};
        }
    SBXCrossover(double rate, double n_c) : crossoverRate(rate), n_c(n_c) {}
    SBXCrossover(const SBXCrossover &) = default;
    SBXCrossover(SBXCrossover &&) = default;
    SBXCrossover &operator=(const SBXCrossover &) = default;
    SBXCrossover &operator=(SBXCrossover &&) = default;
    SBXCrossover() = default;
};
