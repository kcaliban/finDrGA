#include "PrisonersDilemmaGenome.h"

// Use CC & DC from Gen1; CD & DD from Gen2; Init from Gen1
PrisonersDilemmaGenoType PrisonersDilemmaGenome::crossOver(PrisonersDilemmaGenoType gen1, PrisonersDilemmaGenoType gen2, ...) {
  return PrisonersDilemmaGenoType(gen1.strategyCC, gen2.strategyCD,
                                      gen1.strategyDC, gen2.strategyDD,
                                      gen1.strategyInit);
};

PrisonersDilemmaGenoType PrisonersDilemmaGenome::mutate(PrisonersDilemmaGenoType gen, ...) {
  // Change a strategy at random; excluding init from mutation since it does
  // not seem to have a huge effect on performance when many games are simulated
  unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_int_distribution<int> distribution(0, 4);
  int toChange = distribution(generator);
  return PrisonersDilemmaGenoType((toChange == 0) ? not gen.strategyCC : gen.strategyCC,
                                  (toChange == 1) ? not gen.strategyCD : gen.strategyCD,
                                  (toChange == 2) ? not gen.strategyDC : gen.strategyDC,
                                  (toChange == 3) ? not gen.strategyDD : gen.strategyDD,
                                  gen.strategyInit);
};
