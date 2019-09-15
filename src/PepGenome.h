/* Copyright 2019 Fabian Krause
 *
 * PepGA Genome
 *
 * Defines possible mutations and crossover, requires a random engine
*/
#ifndef SRC_PEPGENOME_H_
#define SRC_PEPGENOME_H_
#include <string>
#include <chrono>
#include <random>
#include "lib/Genome.h"
class PepGenome : public Genome<std::string> {
 private:
    std::string alphabet = "ARNDCQEGHILKMFPSTWYV";
    std::mt19937 * mt;

 public:
    PepGenome(std::mt19937 * mt1) {
      mt = mt1;
    }

    std::string crossOver(std::string& , std::string& , ...);
    std::string mutate(std::string&, ...);
};

#endif  // SRC_PEPGENOME_H_
