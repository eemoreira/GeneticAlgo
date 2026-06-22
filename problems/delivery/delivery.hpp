#pragma once

#include <vector>
#include <cmath>
#include <cassert>

struct Delivery {
    double distance_norm, adverse_cond_norm, multiple_deliveries, estimated_time, delivery_person_rating;
    int seg_idx;
    bool is_late;
};


struct DeliveryProblem {
    std::vector<Delivery> deliveries;
    std::vector<double> weights;

    void calcWeights(const std::vector<double>& solution) {
        weights.assign(deliveries.size(), 0.0);
        size_t j = 0;
        for (size_t i = 0; i < 9; i++) {
            double d_max = solution[i * 3 + 0];
            double t_max = solution[i * 3 + 1];
            double f_multi = solution[i * 3 + 2];

            while (j < deliveries.size() && deliveries[j].seg_idx == i) {
                const auto& d = deliveries[j];
                double w_dist = std::exp(-std::max(0.0, d.distance_norm - d_max) * 10);
                double w_adv = std::exp(-std::max(0.0, d.adverse_cond_norm - d_max) * 10);
                double w_multi = 1.0 - (d.multiple_deliveries / 3.0) * (1.0 - f_multi);
                weights[j] = w_dist * w_adv * w_multi;
                j += 1;
            }
        }
        assert(j == deliveries.size());
    }
};

struct DeliveryF1 {
    DeliveryProblem& problem;
    DeliveryF1(DeliveryProblem& problem) : problem(problem) {}
    double operator()(const std::vector<double>& genes) {
        double sum = 0.0;
        double w_sum = 0.0;
        for (size_t i = 0; i < problem.deliveries.size(); i++) {
            sum += problem.weights[i] * problem.deliveries[i].estimated_time;
            w_sum += problem.weights[i];
        }
        return sum / w_sum;
    }
};

struct DeliveryF2 {
    DeliveryProblem& problem;
    DeliveryF2(DeliveryProblem& problem) : problem(problem) {}
    double operator()(const std::vector<double>& genes) {
        double sum = 0.0;
        double w_sum = 0.0;
        for (size_t i = 0; i < problem.deliveries.size(); i++) {
            if (problem.deliveries[i].is_late) {
                sum += problem.weights[i];
            }
            w_sum += problem.weights[i];
        }
        return sum / w_sum;
    }
};

struct DeliveryF3 {
    DeliveryProblem& problem;
    DeliveryF3(DeliveryProblem& problem) : problem(problem) {}
    double operator()(const std::vector<double>& genes) {
        double sum = 0.0;
        double w_sum = 0.0;
        for (size_t i = 0; i < problem.deliveries.size(); i++) {
            sum += problem.weights[i] * problem.deliveries[i].delivery_person_rating;
            w_sum += problem.weights[i];
        }
        return sum / w_sum;
    }
};
