#ifndef STRGN
#define STRGN
#include "../lib/Genome.h"
#include <string>
#include <chrono>
#include <random>
class StringGenome : public Genome<std::string> {
  private:
    std::string alphabet;

  public:
    StringGenome(std::string alphabet1) { alphabet = alphabet1; };

    std::string crossOver(std::string, std::string, ...);

    std::string mutate(std::string, ...);

};

#endif
