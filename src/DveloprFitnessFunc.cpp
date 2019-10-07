/* Copyright 2019 Fabian Krause */
#include "DveloprFitnessFunc.h"

float DveloprFitnessFunc::calculateFitness(std::string & inp, ...) {
  return (-1.0) * poolmgr->getAffinity(inp);
}
