#include <iostream>

#include "generic_algorithm.hpp"
#include "individual.hpp"
#include "crossover.hpp"
#include "mutation.hpp"
#include "rand.hpp"
#include "queens.hpp"
#include <fstream>

signed main(const int argc, const char* argv[]) {

    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <number of queens> <generations> <population_size> <mutation_rate> <crossover_rate> <elitism_selection>" << std::endl;
        return 1;
    }

    int32_t N = std::stoi(argv[1]);
    int32_t generations = std::stoi(argv[2]);
    int32_t populationSize = std::stoi(argv[3]);
    double mutationRate = std::stod(argv[4]);
    double crossoverRate = std::stod(argv[5]);
    uint32_t elitismCharge = std::stoi(argv[6]);
    
    std::cout << "parameters: " << generations << " " << populationSize << " " << mutationRate << " " << crossoverRate << std::endl;

    using GeneType = uint32_t; // permutation

    NQueensFitnessFunction nQueensFitnessFunction;

    auto fitnessFunction = [&](const std::vector<GeneType>& genes) {
        return (double)nQueensFitnessFunction(genes);
    };

    auto generatorFunction = [&]() {
        std::vector<GeneType> permutation(N);
        for (size_t i = 0; i < N; i++) {
            permutation[i] = i;
        }
        std::shuffle(permutation.begin(), permutation.end(), rng);
        return permutation;
    };

    //UniformCrossover crossoverFunction(crossoverRate);
    //OnePointCrossover crossoverFunction(crossoverRate);
    //TwoPointCrossover crossoverFunction(crossoverRate);
    CXCrossover crossoverFunction(crossoverRate);
    //PMXCrossover crossoverFunction(crossoverRate);

    //BitFlipMutation mutationFunction(mutationRate);
    SwapMutation mutationFunction(mutationRate);

    GeneticAlgorithm<GeneType, decltype(fitnessFunction), decltype(generatorFunction), decltype(crossoverFunction), decltype(mutationFunction)> ga(
        populationSize,
        fitnessFunction,
        crossoverFunction,
        mutationFunction,
        generatorFunction,
        elitismCharge
    );

    ga.dimension = N;

    std::ofstream outfile("results.txt", std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "Could not open output file" << std::endl;
        return 1;
    }

   // outfile << "Generations: " << generations << std::endl;
   // outfile << "Population Size: " << populationSize << std::endl;
   // outfile << "Mutation Rate: " << mutationRate << std::endl;
   // outfile << "Crossover Rate: " << crossoverRate << std::endl;
   // outfile << "Elitism Selection: " << elitismCharge << std::endl;

   // outfile << "Generation,BestFitness" << std::endl;

    const int32_t target = N * (N - 1) / 2; // maximum number of good pairs
    int generation = 0;
    double bestFitness = 0.0;
    ga.initPermutationPopulation();
    ga.evaluateFitness();
    while (generation < generations) {
        ga.runSingleGeneration();
        bestFitness = max_element(ga.population.begin(), ga.population.end(), [](const auto& ind1, const auto& ind2) {
            return ind1.fitness < ind2.fitness;
        })->fitness;
        std::cout << "Generation " << generation << ": Best Fitness = " << bestFitness << std::endl;
        outfile << generation << ','  << bestFitness << std::endl;
        generation++;
    }

    std::cout << "Best Fitness: " << bestFitness << std::endl;
    std::cout << "Target Fitness: " << target << std::endl;

    if ((uint32_t)bestFitness == target) {
        std::cout << "Solution found!" << std::endl;

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                if (j == (size_t)ga.population[0].genes[i]) {
                    std::cout << "Q ";
                } else {
                    std::cout << ". ";
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

    } else {
        std::cout << "No solution found." << std::endl;
    }


    outfile << max_element(ga.population.begin(), ga.population.end(), [](const auto& ind1, const auto& ind2) {
        return ind1.fitness < ind2.fitness;
    })->fitness << std::endl;

    return 0;
}


