/* Copyright 2019 Fabian Krause
 *
 * Dvelopr Genome
 *
 * Defines possible mutations and crossover, requires a random engine
*/
#ifndef SRC_DVELOPRGENOME_H_
#define SRC_DVELOPRGENOME_H_
#include <string>
#include <chrono>
#include <random>
#include "lib/Genome.h"
class DveloprGenome : public Genome<std::string> {
 private:
    std::vector<float> probs;
    std::uniform_int_distribution<int> distr;
    std::string alphabet = "ARNDCQEGHILKMFPSTWYV";
    std::mt19937 * mt;

 public:
    DveloprGenome(std::mt19937 * mt1) {
      mt = mt1;
      distr = std::uniform_int_distribution<int>(0, alphabet.size() - 1);
    }

    std::string crossOver(std::string& , std::string& , ...);
    std::string mutate(std::string&, ...);
};

#endif  // SRC_DVELOPRGENOME_H_
