
#pragma once

#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/random.hpp>
#include <vector>
#include "misc.h"

//------------------------------------------------
double runif1(dust::random::xoshiro256plus& rng_state, double min = 0.0, double max = 1.0);

//------------------------------------------------
bool rbernoulli1(dust::random::xoshiro256plus& rng_state, double p);

//------------------------------------------------
int sample2(dust::random::xoshiro256plus& rng_state, int a, int b);

//------------------------------------------------
double rnorm1(dust::random::xoshiro256plus& rng_state, double mean, double sd);

//------------------------------------------------
double rnorm1_interval(dust::random::xoshiro256plus& rng_state, double mean, double sd, double a, double b);

//------------------------------------------------
double rnorm1_pos(dust::random::xoshiro256plus& rng_state, double mean, double sd);
