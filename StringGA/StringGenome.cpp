#include "StringGenome.h"

std::string StringGenome::crossOver(std::string & str1, std::string & str2, ...) {
  // Split at middle (if even, if not then one left to the middle),
  // later could do random point splitting
  int half = str1.size()>>1;
  return str1.substr(0, half) + str2.substr(half, str2.size() - half);
}


std::string StringGenome::mutate(std::string & str1, ...) {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_int_distribution<int> distribution (0, str1.size() - 1);
  std::uniform_int_distribution<int> distributionalphabet (0, alphabet.size() - 1);
  std::string newString = str1;
  newString.at(distribution(generator)) = alphabet.at(distributionalphabet(generator));
  return newString;
}
