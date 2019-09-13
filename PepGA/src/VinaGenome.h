#ifndef VINGN
#define VINGN
#include "lib/Genome.h"
#include <string>
#include <chrono>
#include <random>
class VinaGenome : public Genome<std::string> {
  private:
    std::string alphabet = "ARNDCQEGHILKMFPSTWYV";
    std::mt19937 * mt;

  public:
    VinaGenome(std::mt19937 * mt1) {
      mt = mt1;
    }

    std::string crossOver(std::string& , std::string& , ...);

    std::string mutate(std::string&, ...);

};

#endif
