#ifndef STRFF
#define STRFF
#include "../lib/FitnessFunction.h"
#include <string>
class StringFitnessFunc {
  private:
    std::string goalString;;
  public:
    StringFitnessFunc(std::string goal) { goalString = goal; };
    float calculateFitness(std::string, ...);
};
#endif
