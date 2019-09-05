#ifndef VINGN
#define VINGN
#include "../../lib/Genome.h"
#include <string>
#include <chrono>
#include <random>
class VinaGenome : public Genome<std::string> {
  private:
    std::string alphabet = "ARNDCQEGHILKMFPSTWYV";

  public:
    std::string crossOver(std::string& , std::string& , ...);

    std::string mutate(std::string&, ...);

};

#endif
