#ifndef VINFF
#define VINFF
#include "lib/FitnessFunction.h"
#include "PoolManager/PoolManager.h"
#include <string>
class VinaFitnessFunc {
  private:
    // Pointer to the pool manager
    PoolMGR * poolmgr;

  public:
    VinaFitnessFunc(PoolMGR * poolmgr1) {
      poolmgr = poolmgr1;
    };

    float calculateFitness(std::string &, ...);
};
#endif
