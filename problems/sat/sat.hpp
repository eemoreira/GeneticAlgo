#pragma once

#include <cstdlib>
#include <vector>
#include <cstdint>

struct SatClause {
    std::vector<int32_t> literals;
    bool eval(const std::vector<uint8_t>& assignment) const {
        for (const int32_t& literal : literals) {
            int32_t var_index = abs(literal) - 1;
            bool var_value = assignment[var_index];
            if ((literal > 0 && var_value) || (literal < 0 && !var_value)) {
                return true;
            }
        }
        return false;
    }
};


struct SatFitnessFunction {
    std::vector<SatClause> clauses;

    SatFitnessFunction() = default;
    SatFitnessFunction(const std::vector<SatClause>& cls) : clauses(cls) {}

    uint32_t operator()(const std::vector<uint8_t>& genes) {        uint32_t satisfied_clauses = 0;
        for (const SatClause& clause : clauses) {
            if (clause.eval(genes)) {
                satisfied_clauses += 1;
            }
        }
        return satisfied_clauses;
    }

};

