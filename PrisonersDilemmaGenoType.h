#ifndef PRISGENOTYPE
#define PRISGENOTYPE
class PrisonersDilemmaGenoType {
  public:
    // Strategies for the different cases, 0 ~ C; 1 ~ D
    // Each strategy takes into account own and the others previous action
    bool strategyCC;
    bool strategyCD;
    bool strategyDC;
    bool strategyDD;
    bool strategyInit;

    // Constructor
    PrisonersDilemmaGenoType(bool CC, bool CD, bool DC, bool DD, bool init) {
      strategyCC = CC;
      strategyCD = CD;
      strategyDC = DC;
      strategyDD = DD;
      strategyInit = init;
    };
};

#endif
