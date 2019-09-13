#include "PepFitnessFunc.h"

float PepFitnessFunc::calculateFitness(std::string & inp, ...) {
  return (-1.0) * poolmgr->getAffinity(inp);
}
