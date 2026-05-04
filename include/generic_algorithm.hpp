#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
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

        std::cout << "Evolving population of size " << population.size() << std::endl;

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
