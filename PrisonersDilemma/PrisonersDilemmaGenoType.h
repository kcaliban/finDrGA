#ifndef PRISGENOTYPE
#define PRISGENOTYPE
#include <iostream>
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

    // Destructor
    ~PrisonersDilemmaGenoType() { };

    bool operator == (const PrisonersDilemmaGenoType &Ref) const
    {
        return(this->strategyCC == Ref.strategyCC &&
               this->strategyCD == Ref.strategyCD &&
               this->strategyDC == Ref.strategyDC &&
               this->strategyDD == Ref.strategyDD &&
               this->strategyInit == Ref.strategyInit);
    }

    friend std::ostream& operator<<(std::ostream &strm, const PrisonersDilemmaGenoType &a) {
      return strm << "(Strategy: CC:" << a.strategyCC << ", CD:" << a.strategyCD
                  << ", DC:" << a.strategyDC << ", DD:" << a.strategyDD << "; Init:"
                  << a.strategyInit << ")";
    }
};

#endif
