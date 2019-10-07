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
    std::vector<float> probs;
    std::discrete_distribution<int> distr;
    std::string alphabet = "ARNDCQEGHILKMFPSTWYV";
    std::mt19937 * mt;

 public:
    PepGenome(std::mt19937 * mt1) {
      mt = mt1;
      for (unsigned int i = 0; i < alphabet.size(); i++) {
        if (alphabet.at(i) == 'E' || alphabet.at(i) == 'D') {
          probs.push_back(1.3);
        } else {
          probs.push_back(1.0);
        }
      }
      distr = std::discrete_distribution<int>(probs.begin(), probs.end());
    }

    std::string crossOver(std::string& , std::string& , ...);
    std::string mutate(std::string&, ...);
};

#endif  // SRC_PEPGENOME_H_
