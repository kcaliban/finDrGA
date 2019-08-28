#ifndef PRISGENOM
#define PRISGENOM
#include "Genome.h"
#include "PrisonersDilemmaGenoType.h"
class PrisonersDilemmaGenome : public Genome<PrisonersDilemmaGenoType> {
  public:
    // Use CC & DC from Gen1; CD & DD from Gen2; Init from Gen1
    PrisonersDilemmaGenoType crossOver(PrisonersDilemmaGenoType, PrisonersDilemmaGenoType, ...);

    PrisonersDilemmaGenoType mutate(PrisonersDilemmaGenoType, ...);
};
#endif
