/* Copyright 2019 Fabian Krause */
#ifndef SRC_PEPFITNESSFUNC_H_
#define SRC_PEPFITNESSFUNC_H_
#include <string>
#include "lib/FitnessFunction.h"
#include "PoolManager/PoolManager.h"
class PepFitnessFunc {
 private:
    // Pointer to the pool manager
    PoolMGR * poolmgr;

 public:
    PepFitnessFunc(PoolMGR * poolmgr1) {
      poolmgr = poolmgr1;
    }

    float calculateFitness(std::string &, ...);
};
#endif  // SRC_PEPFITNESSFUNC_H_
