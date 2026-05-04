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
            if (std::uniform_real_distribution(0.0, 1.0)(rng) > crossoverRate) {
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
            if (std::uniform_real_distribution(0.0, 1.0)(rng) > crossoverRate) {
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
            if (std::uniform_real_distribution(0.0, 1.0)(rng) > crossoverRate) {
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
            if (std::uniform_real_distribution(0.0, 1.0)(rng) > crossoverRate) {
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
