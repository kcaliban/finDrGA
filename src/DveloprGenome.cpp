/* Copyright 2019 Fabian Krause */
#include "DveloprGenome.h"

std::string DveloprGenome::crossOver(std::string & str1, std::string & str2, ...) {
  // Excluding direct copy of one of them
  std::uniform_int_distribution<int> distribution(1, str1.size() - 2);
  int split = distribution(*mt);
  return str1.substr(0, split) + str2.substr(split, str2.size() - split);
}


std::string DveloprGenome::mutate(std::string & str1, ...) {
  std::uniform_int_distribution<int> distribution(0, str1.size() - 1);
  std::string newString = str1;
  newString.at(distribution(*mt)) = alphabet.at(distr(*mt));
  return newString;
}
