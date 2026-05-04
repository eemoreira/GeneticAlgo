#include <iostream>

#include "generic_algorithm.hpp"
#include "individual.hpp"
#include "crossover.hpp"
#include "mutation.hpp"
#include "rand.hpp"
#include "sat.hpp"
#include <fstream>
#include <sstream>


struct SatProblem {
  uint32_t num_variables;
  std::vector<SatClause> clauses;

  SatProblem(uint32_t num_vars, const std::vector<SatClause> &cls)
      : num_variables(num_vars), clauses(cls) {}
};

SatProblem read_sat_problem(const std::string& filename) {
    std::vector<SatClause> clauses;
    uint32_t num_variables = 0;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == 'p' || line[0] == '0' || line[0] == '%') {
            continue; // Skip comments and empty lines
        }
        
        std::istringstream iss(line);
        std::vector<int32_t> literals;
        int32_t literal;
        while (iss >> literal) {
            if (literal == 0) {
                break;
            }
            literals.push_back(literal);
            num_variables = std::max(num_variables, static_cast<uint32_t>(std::abs(literal)));
        }
        clauses.push_back({literals});
    }

    return SatProblem(num_variables, clauses);
}


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

    const std::string filename = "problems/sat/sat.cnf";
    SatProblem satProblem = read_sat_problem(filename);

    SatFitnessFunction satFitnessFunction(satProblem.clauses);

    auto fitnessFunction = [&](const std::vector<GeneType>& genes) {
        return satFitnessFunction(genes);
    };

    auto generatorFunction = []() {
        return std::uniform_int_distribution<GeneType>(0, 1)(rng);
    };

    OnePointCrossover crossoverFunction(crossoverRate);

    //BitFlipMutation mutationFunction(mutationRate);
    SwapMutation mutationFunction(mutationRate);

    GeneticAlgorithm<GeneType, decltype(fitnessFunction), decltype(generatorFunction), decltype(crossoverFunction), decltype(mutationFunction)> ga(
        populationSize,
        fitnessFunction,
        crossoverFunction,
        mutationFunction,
        generatorFunction,
        elitismSelection
    );

    ga.dimension = satProblem.num_variables;
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
    while (generation < generations) {
        ga.runSingleGeneration();
        double bestFitness = max_element(ga.population.begin(), ga.population.end(), [](const auto& ind1, const auto& ind2) {
            return ind1.fitness < ind2.fitness;
        })->fitness;
        std::cout << "Generation " << generation << ": Best Fitness = " << bestFitness << std::endl;
        outfile << generation << ','  << bestFitness << std::endl;
        generation++;
    }

    return 0;
}

