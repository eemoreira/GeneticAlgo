#pragma once

#include <vector>

template<typename GeneType>
struct Individual {
    std::vector<GeneType> genes;
    double fitness;

    Individual() : fitness(0.0) {}
    Individual(const std::vector<GeneType>& genes) : genes(genes), fitness(0.0) {}
};

template<typename GeneType>
struct MultiObjectiveIndividual {
    std::vector<GeneType> genes;
    std::vector<double> fitness;
    double crowdingDistance;

    MultiObjectiveIndividual() {}
    MultiObjectiveIndividual(const std::vector<GeneType>& genes) : genes(genes) {}
};
