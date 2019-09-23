#ifndef ser
#define ser
#include <string.h>
#include <vector>
#include <string>
#include <utility>
// Vector of strings
char * serialize(std::vector<std::string> &, unsigned int *);
void deserialize(std::vector<std::string> &, char *, unsigned int);
// Vector of pairs of strings and floats
char * serialize(std::vector<std::pair<std::string, float>> &, unsigned int *);
void deserialize(std::vector<std::pair<std::string, float>> &, char *, unsigned int);
#endif
