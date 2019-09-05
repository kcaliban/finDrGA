#include "VinaFitnessFunc.h"

float VinaFitnessFunc::calculateFitness(std::string & inp, ...) {
  return (-1.0) * poolmgr->getAffinity(inp);
}
