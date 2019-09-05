#include "VinaGenome.h"

std::string VinaGenome::crossOver(std::string & str1, std::string & str2, ...) {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  // Excluding direct copy of one of them
  std::uniform_int_distribution<int> distribution (0, str1.size() - 1);
  int split = distribution(generator);
  return str1.substr(0, split) + str2.substr(split, str2.size() - split);
}


std::string VinaGenome::mutate(std::string & str1, ...) {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_int_distribution<int> distribution (0, str1.size() - 1);
  std::uniform_int_distribution<int> distributionalphabet (0, alphabet.size() - 1);
  std::string newString = str1;
  newString.at(distribution(generator)) = alphabet.at(distributionalphabet(generator));
  return newString;
}
