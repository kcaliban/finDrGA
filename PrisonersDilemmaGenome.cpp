#include "PrisonersDilemmaGenome.h"

// Use CC & DC from Gen1; CD & DD from Gen2; Init from Gen1
PrisonersDilemmaGenoType * PrisonersDilemmaGenome::crossOver(PrisonersDilemmaGenoType gen1, PrisonersDilemmaGenoType gen2, ...) {
  return new PrisonersDilemmaGenoType(gen1.strategyCC, gen2.strategyCD,
                                      gen1.strategyDC, gen2.strategyDD,
                                      gen1.strategyInit);
};

PrisonersDilemmaGenoType * PrisonersDilemmaGenome::mutate(PrisonersDilemmaGenoType gen, ...) {
  // No mutation for now
  return new PrisonersDilemmaGenoType(gen.strategyCC, gen.strategyCD,
                                      gen.strategyDC, gen.strategyDD,
                                      gen.strategyInit);
};
