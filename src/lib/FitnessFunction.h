/* Copyright 2019 iGEM Team Freiburg 2019
 *
 * Defines the fitness function for selection in GA
*/
#ifndef SRC_LIB_FITNESSFUNCTION_H_
#define SRC_LIB_FITNESSFUNCTION_H_
#include "Genome.h"
template <typename GenoType>
class FitnessFunction {
 public:
    virtual float calculateFitness(GenoType, ...) = 0;
};

#endif  // SRC_LIB_FITNESSFUNCTION_H_
