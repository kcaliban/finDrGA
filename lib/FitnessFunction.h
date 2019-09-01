#ifndef FITNS
#define FITNS
#include "Genome.h"
template <typename GenoType>
class FitnessFunction {
  public:
    virtual float calculateFitness(GenoType, ...) = 0;
};

#endif
