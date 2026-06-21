#include <iostream>

#include "generic_algorithm.hpp"
#include "individual.hpp"
#include "crossover.hpp"
#include "mutation.hpp"
#include "rand.hpp"
#include "weighted-queens.hpp"
#include <fstream>

signed main(const int argc, const char* argv[]) {

    if (argc < 8) {
        std::cerr << "Usage: " << argv[0] << " <number of queens> <generations> <population_size> <mutation_rate> <crossover_rate> <elitism_selection> <penalty>" << std::endl;
        return 1;
    }

    int32_t N = std::stoi(argv[1]);
    int32_t generations = std::stoi(argv[2]);
    int32_t populationSize = std::stoi(argv[3]);
    double mutationRate = std::stod(argv[4]);
    double crossoverRate = std::stod(argv[5]);
    uint32_t elitismCharge = std::stoi(argv[6]);
    double penalty = std::stod(argv[7]);
    
    std::cout << "parameters: " << generations << " " << populationSize << " " << mutationRate << " " << crossoverRate << std::endl;

    using GeneType = uint32_t; // permutation

    std::vector<std::vector<double>> a(N, std::vector<double>(N));
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (~i & 1) a[i][j] = sqrt(i * N + j);
            else a[i][j] = log10(i * N + j);
        }
    }

    WeightedQueensFitnessFunction weightedQueensFitnessFunction(penalty, a);

    auto fitnessFunction = [&](const std::vector<GeneType>& genes) {
        return weightedQueensFitnessFunction(genes);
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

    outfile << "Generations: " << generations << std::endl;
    outfile << "Population Size: " << populationSize << std::endl;
    outfile << "Mutation Rate: " << mutationRate << std::endl;
    outfile << "Crossover Rate: " << crossoverRate << std::endl;
    outfile << "Elitism Selection: " << elitismCharge << std::endl;

    outfile << "Generation,BestFitness,AverageFitness" << std::endl;

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
        double averageFitness = std::accumulate(ga.population.begin(), ga.population.end(), 0.0, [](double sum, const auto& ind) {
            return sum + ind.fitness;
        }) / ga.population.size();
        std::cout << "Generation " << generation << ": Best Fitness = " << bestFitness << ", Average Fitness = " << averageFitness << std::endl;
        outfile << generation << ','  << bestFitness << ',' << averageFitness << std::endl;
        generation++;
    }

    std::cout << "Best Fitness: " << bestFitness << std::endl;
    std::cout << "Average Fitness: " << std::accumulate(ga.population.begin(), ga.population.end(), 0.0, [](double sum, const auto& ind) {
        return sum + ind.fitness;
    }) / ga.population.size() << std::endl;

    int bestIndex = std::distance(ga.population.begin(), max_element(ga.population.begin(), ga.population.end(), [](const auto& ind1, const auto& ind2) {
        return ind1.fitness < ind2.fitness;
    }));
    const auto& bestIndividual = ga.population[bestIndex];
    const auto& genes = bestIndividual.genes;
    uint32_t badPairs = 0;
    std::cout << "Best Individual: ";
    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {
            if (abs((int)genes[i] - (int)genes[j]) == abs((int)i - (int)j)) {
                badPairs += 1;
            }
        }
    }
    std::cout << "Number of bad pairs: " << badPairs << std::endl;

    outfile << max_element(ga.population.begin(), ga.population.end(), [](const auto& ind1, const auto& ind2) {
        return ind1.fitness < ind2.fitness;
    })->fitness << std::endl;

    return 0;
}



