#include "PrisonersDilemma.h"
int main(int argc, char *argv[])
{
  GenAlgInst<PrisonersDilemmaGenoType,
             PrisonersDilemmaGenome,
             PrisonersDilemmaFitnessFunc>
    inst;
  // Strategies to train our Genetic Algorithm against
  std::vector<PrisonersDilemmaGenoType> strategies;
  // We want all strategies, i.e. 2^5 different strategies
  for (unsigned int i = 0; i < 32; i++) {
    // Since every strategy is encoded as a boolean, by simply using
    // bitshifting and bit-wise AND we can generate all possible strategies.
    // E.g. the first line:
    // If the bit representation of i has a 0 at its lowest bit (e.g. for all
    // even numbers), we get a cooperate for the CC case, else a deflect.
    strategies.push_back(PrisonersDilemmaGenoType(
            (bool) (i & (1<<0)),
            (bool) (i & (1<<1)),
            (bool) (i & (1<<2)),
            (bool) (i & (1<<3)),
            (bool) (i & (1<<4))));
  }
  std::cout << "Initial generation: [";
  for (auto gen : strategies) {
    std::cout << gen << ",";
  }
  std::cout << "]" << std::endl;
  // Payoff table: number of years "saved" per strategy and players,
  //               i.e. (max number of years) - (actual sentenced years)
  //
  //               Here:
  //               max years = 5
  // Payoff: 1\2| 0(C)  | 1 (D)
  //         ---|----------------
  //       0(C) |   3   |  0
  //       1(D) |   5   |  1
  int payoff[2][2] = {{3, 0}, {5, 1}};
  // Create our fitness function
  PrisonersDilemmaFitnessFunc fitnessFunc =
    PrisonersDilemmaFitnessFunc(strategies, payoff, 100);
  // "Create" a Genome (this is ugly)
  PrisonersDilemmaGenome genome;
  // Run a simulation for 100 rounds, with mutation probability: 0.001
  inst.simulate(genome, fitnessFunc, strategies, 500, 0.001);
}
