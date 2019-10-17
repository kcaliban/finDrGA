/* Copyright 2019 iGEM Team Freiburg 2019
 *
 * finDrGA Genome
 *
 * Defines possible mutations and crossover, requires a random engine
*/
#ifndef SRC_FINDRGAGENOME_H_
#define SRC_FINDRGAGENOME_H_
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include "lib/Genome.h"
class finDrGAGenome : public Genome<std::string> {
 private:
    std::vector<float> probs;
    std::uniform_int_distribution<int> distr;
    std::string alphabet = "ARNDCQEGHILKMFPSTWYV";
    std::mt19937 * mt;

 public:
    finDrGAGenome(std::mt19937 * mt1) {
      mt = mt1;
      distr = std::uniform_int_distribution<int>(0, alphabet.size() - 1);
    }

    std::string crossOver(std::string& , std::string& , ...);
    std::string mutate(std::string&, ...);
};

#endif  // SRC_FINDRGAGENOME_H_
