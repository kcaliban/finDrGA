/* Copyright 2019 Fabian Krause */
#include "finDrGAFitnessFunc.h"

float finDrGAFitnessFunc::calculateFitness(std::string & inp, ...) {
  return (-1.0) * poolmgr->getAffinity(inp);
}
