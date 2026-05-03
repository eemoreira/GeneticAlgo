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

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <generations> <population_size>" << std::endl;
        return 1;
    }

    int32_t generations = atoi(argv[1]);
    int32_t populationSize = atoi(argv[2]);

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

    OnePointCrossover crossoverFunction;

    BitFlipMutation mutationFunction(0.01); // 1% mutation rate

    GeneticAlgorithm<GeneType, decltype(fitnessFunction), decltype(generatorFunction), decltype(crossoverFunction), decltype(mutationFunction)> ga(
        populationSize,
        fitnessFunction,
        crossoverFunction,
        mutationFunction
    );

    ga.dimension = satProblem.num_variables;
    ga.initPopulation();

    // Evaluate initial fitness
    ga.evaluateFitness();

    for (int generation = 0; generation < generations; ++generation) {
        ga.evolve();
        ga.evaluateFitness();
        std::cout << "Generation " << generation << " best fitness: " 
                  << std::max_element(ga.population.begin(), ga.population.end(), [](const auto& a, const auto& b) {
                      return a.fitness < b.fitness;
                  })->fitness << std::endl;
    }

    return 0;
}

