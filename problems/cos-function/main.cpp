#include <cassert>
#include <iostream>
#include "generic_algorithm.hpp"
#include "individual.hpp"
#include "crossover.hpp"
#include "mutation.hpp"
#include "rand.hpp"
#include "cos_function.hpp"
#include <fstream>

signed main(const int argc, const char* argv[]) {
        if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <generations> <population_size> <mutation_rate> <crossover_rate> <elitism_selection>" << std::endl;
        return 1;
    }

    int32_t generations = std::stoi(argv[1]);
    int32_t populationSize = std::stoi(argv[2]);
    double mutationRate = std::stod(argv[3]);
    double crossoverRate = std::stod(argv[4]);
    if (argv[5] != std::string("true") && argv[5] != std::string("false")) {
        std::cerr << "Elitism selection must be 'true' or 'false'" << std::endl;
        return 1;
    }
    bool elitismSelection = argv[5] == std::string("true");
    
    std::cout << "parameters: " << generations << " " << populationSize << " " << mutationRate << " " << crossoverRate << std::endl;

    using GeneType = uint8_t; // Binary representation (not using bool because std::bool sucks)

    CosFitnessFunction cosFitnessFunction;
    auto DecodeBinaryElement = [&](const std::vector<GeneType>& genes) {
        double value = 0.0;
        for (size_t i = 0; i < genes.size(); ++i) {
            value += genes[i] * (1 << i);
        }
        // real value in range [-2, 2]
        value = (value / (double)((1 << genes.size()) - 1)) * 4 - 2;
        return value;
    };

    auto fitnessFunction = [&](const std::vector<GeneType>& genes) {
        double x = DecodeBinaryElement(genes);
        return std::max(0.0, cosFitnessFunction(x));
    };

    auto generatorFunction = []() {
        return std::uniform_int_distribution<GeneType>(0, 1)(rng);
    };

    OnePointCrossover crossoverFunction(crossoverRate);

    BitFlipMutation mutationFunction(mutationRate);
    //SwapMutation mutationFunction(mutationRate);

    GeneticAlgorithm<GeneType, decltype(fitnessFunction), decltype(generatorFunction), decltype(crossoverFunction), decltype(mutationFunction)> ga(
        populationSize,
        fitnessFunction,
        crossoverFunction,
        mutationFunction,
        generatorFunction,
        elitismSelection
    );

    ga.dimension = 30;
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

    return 0;
}
