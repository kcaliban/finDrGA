/* Copyright 2019 Fabian Krause
 *
 * PepGA Fitness function
 *
 * Returns negated binding affinities for the use of the genetic algorithm,
 * requires a PoolMGR for the affinities
*/
#ifndef SRC_PEPFITNESSFUNC_H_
#define SRC_PEPFITNESSFUNC_H_
#include <string>
#include "lib/FitnessFunction.h"
#include "PoolManager/PoolManager.h"
class PepFitnessFunc {
 private:
    PoolMGR * poolmgr;

 public:
    PepFitnessFunc(PoolMGR * poolmgr1) {
      poolmgr = poolmgr1;
    }

    float calculateFitness(std::string &, ...);
};

#endif  // SRC_PEPFITNESSFUNC_H_
