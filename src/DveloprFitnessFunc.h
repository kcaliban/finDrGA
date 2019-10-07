/* Copyright 2019 Fabian Krause
 *
 * Dvelopr Fitness function
 *
 * Returns negated binding affinities for the use of the genetic algorithm,
 * requires a PoolMGR for the affinities
*/
#ifndef SRC_DVELOPRFITNESSFUNC_H_
#define SRC_DVELOPRFITNESSFUNC_H_
#include <string>
#include "lib/FitnessFunction.h"
#include "PoolManager/PoolManager.h"
class DveloprFitnessFunc {
 private:
    PoolMGR * poolmgr;

 public:
    DveloprFitnessFunc(PoolMGR * poolmgr1) {
      poolmgr = poolmgr1;
    }

    float calculateFitness(std::string &, ...);
};

#endif  // SRC_DVELOPRFITNESSFUNC_H_
