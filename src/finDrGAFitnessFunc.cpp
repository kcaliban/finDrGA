/* Copyright 2019 iGEM Team Freiburg 2019 */
#include "finDrGAFitnessFunc.h"

float finDrGAFitnessFunc::calculateFitness(std::string & inp, ...) {
  return (-1.0) * poolmgr->getAffinity(inp);
}
