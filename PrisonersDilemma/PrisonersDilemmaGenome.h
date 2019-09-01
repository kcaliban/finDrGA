#ifndef PRISGENOM
#define PRISGENOM
#include "../lib/Genome.h"
#include "PrisonersDilemmaGenoType.h"
#include <random>
#include <chrono>
class PrisonersDilemmaGenome : public Genome<PrisonersDilemmaGenoType> {
  public:
    // Use CC & DC from Gen1; CD & DD from Gen2; Init from Gen1
    PrisonersDilemmaGenoType crossOver(PrisonersDilemmaGenoType, PrisonersDilemmaGenoType, ...);

    PrisonersDilemmaGenoType mutate(PrisonersDilemmaGenoType, ...);
};
#endif
