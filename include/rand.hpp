#pragma once

#include <random>
#include <chrono>

static std::mt19937 rng((uint32_t)std::chrono::steady_clock::now().time_since_epoch().count());
