/* Copyright 2019 iGEM Team Freiburg 2019
 *
 * String and (string, float) pair serialization required for communication
 * between master and workers
 *
*/
#ifndef SRC_SERIALIZATION_SERIALIZATION_H_
#define SRC_SERIALIZATION_SERIALIZATION_H_
#include <string.h>
#include <vector>
#include <string>
#include <utility>
// Vector of strings
char * serialize(std::vector<std::string> &, unsigned int *);
void deserialize(std::vector<std::string> &, char *, unsigned int);
// Vector of pairs of strings and floats
char * serialize(std::vector<std::pair<std::string, float>> &, unsigned int *);
void deserialize(std::vector<std::pair<std::string, float>> &,
                 char *, unsigned int);
#endif  // SRC_SERIALIZATION_SERIALIZATION_H_
