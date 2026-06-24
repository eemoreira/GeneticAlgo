#include "generic_algorithm.hpp"
#include "crossover.hpp"
#include "mutation.hpp"
#include "rand.hpp"
#include "delivery.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <stdexcept>
#include <functional>


static std::vector<std::string> split(const std::string& line, char delim = ',') {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string item;

    while (std::getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

static std::unordered_map<std::string, int> build_header_map(const std::string& header_line) {
    auto cols = split(header_line);
    std::unordered_map<std::string, int> mp;

    for (int i = 0; i < (int)cols.size(); i++) {
        std::string key = cols[i];
        key.erase(0, key.find_first_not_of(" \t\r\n"));
        key.erase(key.find_last_not_of(" \t\r\n") + 1);

        mp[key] = i;
    }

    return mp;
}

DeliveryProblem load_from_csv(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Nao foi possivel abrir o CSV");
    }

    std::string line;

    std::getline(file, line);
    auto idx = build_header_map(line);

    auto require = [&](const std::string& col) {
        if (!idx.count(col)) {
            throw std::runtime_error("Coluna ausente no CSV: " + col);
        }
    };

    require("dist_norm");
    require("adverso_norm");
    require("multiple_deliveries");
    require("seg_idx");
    require("tempo_estimado");
    require("Delivery_person_Ratings");
    require("is_late");

    DeliveryProblem problem;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        auto cols = split(line);

        Delivery d;

        d.distance_norm = std::stod(cols[idx["dist_norm"]]);
        d.adverse_cond_norm = std::stod(cols[idx["adverso_norm"]]);
        d.multiple_deliveries = std::stod(cols[idx["multiple_deliveries"]]);
        d.seg_idx = std::stoi(cols[idx["seg_idx"]]);
        d.estimated_time = std::stod(cols[idx["tempo_estimado"]]);
        d.delivery_person_rating = std::stod(cols[idx["Delivery_person_Ratings"]]);
        d.is_late = std::stoi(cols[idx["is_late"]]) != 0;

        problem.deliveries.push_back(d);
    }

    sort(problem.deliveries.begin(), problem.deliveries.end(), [](const Delivery& a, const Delivery& b) {
        return a.seg_idx < b.seg_idx;
    });
    
    return problem;
}

signed main(int argc, const char* argv[]) {

    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <generations> <population_size> <mutation_rate> <crossover_rate>" << std::endl;
        return 1;
    }

    const std::string filename = "problems/delivery/delivery_dataset.csv";
    int32_t generations = std::stoi(argv[1]);
    int32_t populationSize = std::stoi(argv[2]);
    double mutationRate = std::stod(argv[3]);
    double crossoverRate = std::stod(argv[4]);

    DeliveryProblem problem;
    try {
        problem = load_from_csv(filename);
    } catch (const std::exception& ex) {
        std::cerr << "Erro: " << ex.what() << std::endl;
        return 1;
    }

    DeliveryF1 f1(problem);
    DeliveryF2 f2(problem);
    DeliveryF3 f3(problem);

    using FitnessFun = std::function<double(const std::vector<double>&)>;
    using ComparisonFun = std::function<bool(double, double)>;

    FitnessFun F1 = [&](const std::vector<double>& genes) { return f1(genes); };
    FitnessFun F2 = [&](const std::vector<double>& genes) { return f2(genes); };
    FitnessFun F3 = [&](const std::vector<double>& genes) { return f3(genes); };
    std::vector<FitnessFun> fitnessFunctions = {F1, F2, F3};
    ComparisonFun cmp1 = [&](double a, double b) { return a < b; }; // minimizar tempo
    ComparisonFun cmp2 = [&](double a, double b) { return a < b; }; // minimizar atraso
    ComparisonFun cmp3 = [&](double a, double b) { return a > b; }; // maximizar avaliação
    std::vector<ComparisonFun> comparisonFunctions = {cmp1, cmp2, cmp3};


    using GeneType = double; // [0, 1] para cada gene
    
    auto generatorFunction = [&]() {
        return std::uniform_real_distribution<>(0.0, 1.0)(rng);
    };

    SBXCrossover crossoverFunction(crossoverRate, 10.0);
    PolynomialMutation mutationFunction(mutationRate, 10.0);

    NSGA2<
        GeneType,
        decltype(generatorFunction),
        decltype(crossoverFunction),
        decltype(mutationFunction),
        FitnessFun,
        ComparisonFun
    > ga(
        populationSize,
        crossoverFunction,
        mutationFunction,
        generatorFunction,
        fitnessFunctions,
        comparisonFunctions
    );

    ga.dimension = 27; // 9 segmentos * 3 parâmetros cada

    ga.initPopulation();
    while (generations--) {
        std::cout << "Generations left = " << generations + 1 << " \r" << std::flush;
        for (auto& ind : ga.population) {
            problem.calcWeights(ind.genes);
            ind.fitness[0] = F1(ind.genes);
            ind.fitness[1] = F2(ind.genes);
            ind.fitness[2] = F3(ind.genes);
        }
        ga.evolve();
    }

    std::ofstream outfile("problems/delivery/output.csv");
    outfile << "f1,f2,f3\n";

    auto pareto = ga.getFronts()[0];
    std::cout << "Pareto front size: " << pareto.size() << std::endl;

    for (auto& ind : pareto) {
        problem.calcWeights(ind.genes);
        ind.fitness[0] = F1(ind.genes);
        ind.fitness[1] = F2(ind.genes);
        ind.fitness[2] = F3(ind.genes);
        outfile << ind.fitness[0] << "," << ind.fitness[1] << "," << ind.fitness[2] << "\n";
    }

    return 0;
}
