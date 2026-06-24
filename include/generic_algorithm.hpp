#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include "individual.hpp"
#include "rand.hpp"

template<
    typename GeneType,
    typename FitnessFunction,
    typename GeneratorFunction,
    typename CrossoverFunction,
    typename MutationFunction
>
struct GeneticAlgorithm {
    using IndividualType = Individual<GeneType>;

    std::vector<IndividualType> population;
    int32_t dimension;
    FitnessFunction fitnessFunction;
    CrossoverFunction crossoverFunction;
    MutationFunction mutationFunction;
    GeneratorFunction generatorFunction;
    uint32_t elitismCharge;

    GeneticAlgorithm(
        const uint32_t& populationSize,
        FitnessFunction fitnessFunc,
        CrossoverFunction crossoverFunc,
        MutationFunction mutationFunc,
        GeneratorFunction generatorFunc,
        uint32_t elitism
    ) : population(populationSize),
        fitnessFunction(fitnessFunc),
        crossoverFunction(crossoverFunc),
        mutationFunction(mutationFunc),
        generatorFunction(generatorFunc),
        elitismCharge(elitism) {}

    void initPopulation() {
        for (auto& ind : population) {
            ind.genes.resize(dimension);
            randomize(ind);
        }
    }

    void initPermutationPopulation() {
        for (auto& ind : population) {
            ind.genes.resize(dimension);
            randomizePermutation(ind);
        }
    }

    void randomizePermutation(IndividualType& individual) {
        individual.genes = generatorFunction();
    }

    void randomize(IndividualType& individual) {
        for (auto& g : individual.genes) {
            g = generatorFunction();
        }
    }

    void evaluateFitness() {
        for (auto& ind : population) {
            ind.fitness = (double) fitnessFunction(ind.genes);
        }
    }

    void run(int generations) {
        for (int i = 0; i < generations; ++i) {
            evolve();
            evaluateFitness();
        }
    }

    void runSingleGeneration() {
        evolve();
        evaluateFitness();
    }

    size_t selectParent() {
        double totalFitness = 0.0;
        for (const auto& ind : population) {
            totalFitness += ind.fitness;
        }

        double r = std::uniform_real_distribution<>(0.0, totalFitness)(rng);
        double now = 0.0;
        for (size_t i = 0; i < population.size(); ++i) {
            const auto& ind = population[i];
            now += ind.fitness;
            if (now >= r) {
                return i;
            }
        }
        return population.size() - 1; // should not happen if totalFitness > 0
    }

    void evolve() {
        std::vector<IndividualType> newPopulation;
        newPopulation.reserve(population.size());

        std::nth_element(population.begin(), population.begin() + elitismCharge, population.end(), [](const auto& a, const auto& b) {
                return a.fitness > b.fitness;
        });
        for (uint32_t i = 0; i < elitismCharge; i++) {
            newPopulation.emplace_back(population[i]);
        }


        while (newPopulation.size() < population.size()) {
            size_t p1 = selectParent();
            size_t p2 = selectParent();

            while (p2 == p1) {
                p2 = selectParent();
            }

            const auto& parent1 = population[p1];
            const auto& parent2 = population[p2];
            auto [c1, c2] = crossoverFunction(parent1, parent2);

            mutationFunction(c1);
            newPopulation.push_back(std::move(c1));

            if (newPopulation.size() < population.size()) {
                mutationFunction(c2);
                newPopulation.push_back(std::move(c2));
            }
        }

        swap(population, newPopulation);
    }

};

template<
    typename GeneType,
    typename GeneratorFunction,
    typename CrossoverFunction,
    typename MutationFunction,
    typename FitnessFunction,
    typename DominatesOnFunction
>
struct NSGA2 {
    using IndividualType = MultiObjectiveIndividual<GeneType>;

    std::vector<IndividualType> population;
    std::vector<FitnessFunction> fitnessFunctions;
    int32_t dimension;
    size_t populationSize;
    CrossoverFunction crossoverFunction;
    MutationFunction mutationFunction;
    GeneratorFunction generatorFunction;
    std::vector<DominatesOnFunction> dominatesOnFunction;

    NSGA2(
        const size_t& popSize,
        CrossoverFunction crossoverFunc,
        MutationFunction mutationFunc,
        GeneratorFunction generatorFunc,
        const std::vector<FitnessFunction>& fitnessFuncs,
        const std::vector<DominatesOnFunction>& dominatedOnFunc
    ) : population(popSize),
        populationSize(popSize),
        crossoverFunction(crossoverFunc),
        mutationFunction(mutationFunc),
        generatorFunction(generatorFunc),
        dominatesOnFunction(dominatedOnFunc),
        fitnessFunctions(fitnessFuncs) {}

    void initPopulation() {
        for (auto& ind : population) {
            ind.genes.resize(dimension);
            ind.fitness.resize(fitnessFunctions.size());
            randomize(ind);
        }
    }

    void randomize(IndividualType& individual) {
        for (auto& g : individual.genes) {
            g = generatorFunction();
        }
    }

    bool dominates(const IndividualType& a, const IndividualType& b) {
        bool betterInAtLeastOne = false;
        assert(a.fitness.size() == b.fitness.size());
        for (size_t i = 0; i < a.fitness.size(); ++i) {
            assert(i < dominatesOnFunction.size());
            const auto& cmp = dominatesOnFunction[i];

            if (cmp(b.fitness[i], a.fitness[i])) {
                return false;
            }

            if (cmp(a.fitness[i], b.fitness[i])) {
                betterInAtLeastOne = true;
            }
        }
        return betterInAtLeastOne;
    }

    size_t selectParent() {
        int idx1 = std::uniform_int_distribution<>(0, population.size() - 1)(rng);
        int idx2 = std::uniform_int_distribution<>(0, population.size() - 1)(rng);

        if (dominates(population[idx1], population[idx2])) {
            return idx1;
        } else if (dominates(population[idx2], population[idx1])) {
            return idx2;
        } else {
            return population[idx1].crowdingDistance > population[idx2].crowdingDistance ? idx1 : idx2;
        }
    }

    std::vector<IndividualType> getNewPopulation() {
        std::vector<IndividualType> newPopulation;
        newPopulation.reserve(populationSize);

        while (newPopulation.size() < populationSize) {
            size_t p1 = selectParent();
            size_t p2 = selectParent();

            while (p2 == p1) {
                p2 = selectParent();
            }

            const auto& parent1 = population[p1];
            const auto& parent2 = population[p2];
            auto [c1, c2] = crossoverFunction(parent1, parent2);

            mutationFunction(c1);
            newPopulation.push_back(std::move(c1));

            if (newPopulation.size() < populationSize) {
                mutationFunction(c2);
                newPopulation.push_back(std::move(c2));
            }
        }
        return newPopulation;
    }

    void assignCrowdingDistances(std::vector<IndividualType>& front) {
        for (auto& ind : front) {
            ind.crowdingDistance = 0.0;
        }

        for (size_t m = 0; m < front[0].fitness.size(); ++m) {
            std::sort(front.begin(), front.end(), [m](const auto& a, const auto& b) {
                return a.fitness[m] < b.fitness[m];
            });

            front[0].crowdingDistance = front.back().crowdingDistance = std::numeric_limits<double>::infinity();

            double fMin = front[0].fitness[m];
            double fMax = front.back().fitness[m];

            if (fMax - fMin == 0) continue;

            for (size_t i = 1; i + 1 < front.size(); ++i) {
                front[i].crowdingDistance += (front[i + 1].fitness[m] - front[i - 1].fitness[m]) / (fMax - fMin);
            }
        }
    }

    std::vector<std::vector<IndividualType>> getFronts() {
        std::vector<std::vector<size_t>> fronts;
        std::vector<std::vector<int>> dom(population.size());
        std::vector<int> domCnt(population.size(), 0);

        for (size_t i = 0; i < population.size(); i++) {
            for (size_t j = 0; j < population.size(); j++) {
                if (i != j && dominates(population[i], population[j])) {
                    dom[i].emplace_back(j);
                    domCnt[j]++;
                }
            }
        }

        for (size_t i = 0; i < population.size(); i++) {
            if (domCnt[i] == 0) {
                if (fronts.empty()) {
                    fronts.emplace_back();
                }
                fronts[0].emplace_back(i);
            }
        }

        std::vector<std::vector<IndividualType>> resultFronts;
        for (const auto& front : fronts) {
            std::vector<IndividualType> frontIndividuals;
            for (size_t idx : front) {
                frontIndividuals.emplace_back(population[idx]);
            }
            resultFronts.emplace_back(std::move(frontIndividuals));
        }

        while (true) {
            std::vector<size_t> nextFront;
            for (const auto& i : fronts.back()) {
                for (size_t j : dom[i]) {
                    if (--domCnt[j] == 0) {
                        nextFront.emplace_back(j);
                    }
                }
            }
            if (nextFront.empty()) {
                break;
            }
            fronts.emplace_back(std::move(nextFront));
            std::vector<IndividualType> frontIndividuals;
            for (size_t idx : fronts.back()) {
                frontIndividuals.emplace_back(population[idx]);
            }
            resultFronts.emplace_back(std::move(frontIndividuals));
        }

        return resultFronts;
    }

    void evolve() {
        assignCrowdingDistances(population);
        auto offspring = getNewPopulation();
        population.insert(population.end(), offspring.begin(), offspring.end());

        auto fronts = getFronts();
        std::vector<IndividualType> nxtPopulation;

        for (auto& front : fronts) {
            if (nxtPopulation.size() + front.size() <= populationSize) {
                nxtPopulation.insert(nxtPopulation.end(), front.begin(), front.end());
            } else {
                assignCrowdingDistances(front);
                std::sort(front.begin(), front.end(), [](const auto& a, const auto& b) {
                    return a.crowdingDistance > b.crowdingDistance;
                });
                for (size_t i = 0; nxtPopulation.size() < populationSize; ++i) {
                    assert(i < front.size());
                    nxtPopulation.emplace_back(front[i]);
                }
                break;
            }
        }
        swap(population, nxtPopulation);
    }
};

