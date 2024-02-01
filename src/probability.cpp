
#include "probability.h"


//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(dust::random::xoshiro256plus& rng_state, double min, double max) {
  return dust::random::uniform(rng_state, min, max);
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal
// probability
int sample2(dust::random::xoshiro256plus& rng_state, int a, int b) {
  return floor(runif1(rng_state, a, b + 1));
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(dust::random::xoshiro256plus& rng_state, double mean, double sd) {
  return dust::random::normal(rng_state, mean, sd);
}

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(dust::random::xoshiro256plus& rng_state, double mean, double sd, double a, double b) {
  
  // draw raw value relative to a
  double ret = rnorm1(rng_state, mean, sd) - a;
  
  // reflect off boundaries at 0 and (b - a)
  if (ret < 0 || ret > (b - a)) {
    
    // use multiple reflections to bring into range [-(b - a), 2(b - a)]
    if (ret < -(b - a)) {
      int n_double_intervals = floor(-ret / (b - a)) / 2;
      ret += 2 * (b - a) * (n_double_intervals + 1);
    } else if (ret > 2*(b - a)) {
      int n_double_intervals = floor(ret / (b - a) - 1) / 2;
      ret -= 2 * (b - a) * (n_double_intervals + 1);
    }
    
    // use one more reflection to bring into range [0, (b - a)]
    if (ret < 0) {
      ret = -ret;
    }
    if (ret > (b - a)) {
      ret = 2*(b - a) - ret;
    }
  }
  
  // no longer relative to a
  ret += a;
  
  // don't let ret equal exactly a or b
  if (ret == a) {
    ret += UNDERFLO_DOUBLE;
  } else if (ret == b) {
    ret -= UNDERFLO_DOUBLE;
  }
  
  return ret;
}

//------------------------------------------------
// draw from reflected normal distribution about 0
double rnorm1_pos(dust::random::xoshiro256plus& rng_state, double mean, double sd) {
  
  // draw raw value
  double ret = rnorm1(rng_state, mean, sd);
  
  // reflect
  if (ret < 0) {
    ret = -ret;
  }
  
  // don't let ret equal exactly 0
  if (ret == 0) {
    ret += UNDERFLO_DOUBLE;
  }
  
  return ret;
}
