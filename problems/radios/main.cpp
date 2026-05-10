#include "generic_algorithm.hpp"
#include "crossover.hpp"
#include "mutation.hpp"
#include "rand.hpp"
#include "radios.h"
#include <fstream>

signed main(const int argc, const char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <generations> <population_size> <mutation_rate> <crossover_rate> <elitism_selection> <lambda>" << std::endl;
        return 1;
    }

    int32_t generations = std::stoi(argv[1]);
    int32_t populationSize = std::stoi(argv[2]);
    double mutationRate = std::stod(argv[3]);
    double crossoverRate = std::stod(argv[4]);
    int32_t elitismSelection = std::stoi(argv[5]);
    double lambda = std::stod(argv[6]);
    
    std::cout << "parameters: " << generations << " " << populationSize << " " << mutationRate << " " << crossoverRate << std::endl;

    using GeneType = uint8_t; // Binary representation (not using bool because std::bool sucks)

    RadioFitnessFunction radioFitnessFunction(lambda);

    auto generatorFunction = []() {
        return std::uniform_int_distribution<GeneType>(0, 1)(rng);
    };

    OnePointCrossover crossoverFunction(crossoverRate);
    BitFlipMutation mutationFunction(mutationRate);
    
    GeneticAlgorithm<GeneType, decltype(radioFitnessFunction), decltype(generatorFunction), decltype(crossoverFunction), decltype(mutationFunction)> ga(
        populationSize,
        radioFitnessFunction,
        crossoverFunction,
        mutationFunction,
        generatorFunction,
        elitismSelection
    );

    ga.dimension = 10; // (5 for standard, 5 for luxury)
    ga.initPopulation();

    int generation = 0;
    // write the best fitness of each generation to a file
    std::ofstream outfile("output.txt");
    if (!outfile.is_open()) {
        std::cerr << "Could not open output file" << std::endl;
        return 1;
    }

    outfile << "Generations: " << generations << std::endl;
    outfile << "Population Size: " << populationSize << std::endl;
    outfile << "Mutation Rate: " << mutationRate << std::endl;
    outfile << "Crossover Rate: " << crossoverRate << std::endl;
    outfile << "Elitism Selection: " << (elitismSelection ? "true" : "false") << std::endl;

    outfile << "Generation,BestFitness" << std::endl;
    ga.evaluateFitness();
    double maxFitness = 0.0;
    while (generation < generations) {
        ga.runSingleGeneration();
        double bestFitness = max_element(ga.population.begin(), ga.population.end(), [](const auto& ind1, const auto& ind2) {
            return ind1.fitness < ind2.fitness;
        })->fitness;
        maxFitness = std::max(maxFitness, bestFitness);
        std::cout << "Generation " << generation << ": Best Fitness = " << bestFitness << std::endl;
        outfile << generation << ','  << bestFitness << std::endl;
        generation++;
    }

    std::cout << "Max Fitness: " << maxFitness << std::endl;
    std::cout << "Fraction from optimal (with no constraints): " << maxFitness / (24 * 30 + 16 * 40) << std::endl;
    auto [standard, luxury] = radioFitnessFunction.decode(ga.population[0].genes);
    std::cout << "Best solution: standard = " << standard << ", luxury = " << luxury << std::endl;
    std::cout << "Penalty: " << radioFitnessFunction.penalty(std::max(0, standard - 24), std::max(0, luxury - 16), std::max(0, standard + 2*luxury - 40)) << std::endl;

    return 0;
}
