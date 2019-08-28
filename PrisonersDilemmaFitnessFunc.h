#ifndef PRISFITNESS
#define PRISFITNESS
#include "PrisonersDilemmaGenoType.h"
#include "FitnessFunction.h"
#include <vector>
#include <string.h>
class PrisonersDilemmaFitnessFunc : public FitnessFunction<PrisonersDilemmaGenoType>
{
  public:
    PrisonersDilemmaFitnessFunc(
        std::vector<PrisonersDilemmaGenoType> strategies1,
        int payoff1[2][2], int n1) {
      strategies = strategies1;
      n = n1;
      memcpy(payoff, payoff1, sizeof(payoff));
    };

    // using FitnessFunction<PrisonersDilemmaGenoType>::calculateFitness;
    float calculateFitness(PrisonersDilemmaGenoType, ...);
  private:
    std::vector<PrisonersDilemmaGenoType> strategies;
    int payoff[2][2];
    int n;
    bool * simGame(PrisonersDilemmaGenoType, PrisonersDilemmaGenoType, bool[2]);
};

#endif
