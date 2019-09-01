#include "StringFitnessFunc.h"
#include <iostream>

float StringFitnessFunc::calculateFitness(std::string input, ...) {
  // Naive Hamming distance
  unsigned int score = input.size();
  for (unsigned int i = 0; i < input.size(); i++) {
    if (input.at(i) != goalString.at(i)) {
      score--;
    }
  }
  if (score == input.size()) {
    std::cout << "\tGoal String found!" << std::endl;
    exit(0);
  }

  return (float) score;
}
