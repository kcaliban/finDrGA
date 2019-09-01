#include "PrisonersDilemmaFitnessFunc.h"

void PrisonersDilemmaFitnessFunc::simGame(PrisonersDilemmaGenoType strat1, PrisonersDilemmaGenoType strat2, bool prevGame[2], bool * newStrat) {
  // Given a previous game and two strategies, return the next moves
  // bool * newStrat = new bool[2];
  if (prevGame[0] == false) {
    if (prevGame[1] == false) {
      newStrat[0] = strat1.strategyCC;
      newStrat[1] = strat2.strategyCC;
    } else {
      newStrat[0] = strat1.strategyCD;
      newStrat[1] = strat2.strategyCD;
    }
  } else {
    if (prevGame[1] == false) {
      newStrat[0] = strat1.strategyDC;
      newStrat[1] = strat2.strategyDC;
    } else {
      newStrat[0] = strat1.strategyDD;
      newStrat[1] = strat2.strategyDD;
    }
  }
}

float PrisonersDilemmaFitnessFunc::calculateFitness(PrisonersDilemmaGenoType gen, ...) {
  // Calculate Fitness for Gen by simulating n games against each strategy
  // Payoff: 1\2| 0(C)  | 1 (D)
  //         ---|----------------
  //       0(C) |       |
  //       1(D) |       |
  // Payoff is given just for gen
  float fitness = 0;
  for (auto strat : strategies) {
    bool prevGame[2] = {gen.strategyInit, strat.strategyInit};
    for (int i = 0; i < n; i++) {
      bool game[2];
      simGame(gen, strat, prevGame, game);
      if (game[0] == false) {
        if (game[1] == false) {
          fitness += payoff[0][0];
        } else {
          fitness += payoff[0][1];
        }
      } else {
        if (game[1] == false) {
          fitness += payoff[1][0];
        } else {
          fitness += payoff[1][1];
        }
      }
    }
  }

  return (float)fitness / (float)n;
}
