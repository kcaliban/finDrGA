#include "Serialization.h"

using namespace std;

char * serialize(vector<string> &v, unsigned int *size) {
  // STR1\0STR2\0...STRN\0
  unsigned int totalSize = 0;

  // Length of strings, +1 for \0
  for (auto str : v) {
    totalSize += str.size() + 1;
  }

  char * buffer = new char[totalSize];
  unsigned int bufpt = 0;

  for (auto str : v) {
    for (unsigned int j = 0; j < str.size(); j++) {
      buffer[bufpt++] = str[j];
    }
    buffer[bufpt++] = '\0';
  }

  * size = totalSize;

  return buffer;
}

void deserialize(vector<string> &restore,  char* buffer, unsigned int size) {
  string curStr;
  for (unsigned int j = 0; j < size; j++) {
    const char c = buffer[j];
    if (c != '\0') {
      curStr.push_back(c);
    } else {
      restore.push_back(curStr);
      curStr.clear();
    }
  }
}

/********************************************/
char * serialize(vector<pair<string, float>> &v, unsigned int *size) {
  // Floats can have \0, \3 (fittingly "end of text") is unused
  // STR1\3FL1\3STR2\3FL2\3...STRN\3FLN\3
  unsigned int totalSize = 0;

  // Length of string, +1 for \0, size of float, +1 for \0
  for (auto p : v) {
    totalSize += p.first.size() + 1 + sizeof(float) + 1;
  }

  char * buffer = new char[totalSize];
  unsigned int bufpt = 0;

  for (auto p : v) {
    for (unsigned int j = 0; j < p.first.size(); j++) {
      buffer[bufpt++] = p.first[j];
    }
    buffer[bufpt++] = '\3';

    char flt[sizeof(float)];
    float val = p.second;
    memcpy(&flt, &val, sizeof(float));
    for (unsigned int j = 0; j < sizeof(float); j++) {
      buffer[bufpt++] = flt[j];
    }
    buffer[bufpt++] = '\3';
  }

  * size = totalSize;
  return buffer;
}

void deserialize(vector<pair<string, float>> &restore, char * buffer, unsigned int size) {
  bool firstterm = false;
  string curStr;
  char flt[sizeof(float)];
  unsigned int fltptr = 0;

  for (unsigned int j = 0; j < size; j++) {
    const char c = buffer[j];
    if (c != '\3') {
      if (!firstterm) {
        curStr.push_back(c);
      }   // Read in string
      else {
        flt[fltptr++] = buffer[j];
      }   // Read in float
    } else {
      if (firstterm) {
        float number;
        memcpy(&number, &flt, sizeof(float));
        restore.push_back(make_pair(curStr, number));

        curStr.clear();
        fltptr = 0;
        firstterm = false;
      } else {
        firstterm = true;
      }
    }
  }
}
