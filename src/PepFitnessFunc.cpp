/* Copyright 2019 Fabian Krause */
#include "PepFitnessFunc.h"

float PepFitnessFunc::calculateFitness(std::string & inp, ...) {
  return (-1.0) * poolmgr->getAffinity(inp);
}
