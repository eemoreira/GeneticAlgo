#pragma once

#include <vector>

template<typename GeneType>
struct Individual {
    std::vector<GeneType> genes;
    double fitness;

    Individual() : fitness(0.0) {}
    Individual(const std::vector<GeneType>& genes) : genes(genes), fitness(0.0) {}
};
