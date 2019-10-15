/* Copyright 2019 Fabian Krause
 *
 * finDrGA Fitness function
 *
 * Returns negated binding affinities for the use of the genetic algorithm,
 * requires a PoolMGR for the affinities
*/
#ifndef SRC_FINDRGAFITNESSFUNC_H_
#define SRC_FINDRGAFITNESSFUNC_H_
#include <string>
#include "lib/FitnessFunction.h"
#include "PoolManager/PoolManager.h"
class finDrGAFitnessFunc {
 private:
    PoolMGR * poolmgr;

 public:
    finDrGAFitnessFunc(PoolMGR * poolmgr1) {
      poolmgr = poolmgr1;
    }

    float calculateFitness(std::string &, ...);
};

#endif  // SRC_FINDRGAFITNESSFUNC_H_
