// Copyright 2019 iGEM Team Freiburg 2019
#include "Serialization.h"
#include <iostream>

using std::string;
using std::vector;
using std::pair;

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
  // STR1\0FL1STR2\0FL2...STRN\0FLN
  unsigned int totalSize = 0;

  // Length of string, +1 for \0, size of float
  for (auto it = v.begin(); it != v.end(); it++) {
    totalSize += it->first.size() + 1 + sizeof(float);
  }
  totalSize++;

  char * buffer = new char[totalSize];
  unsigned int bufpt = 0;

  for (auto it = v.begin(); it != v.end(); it++) {
    for (unsigned int j = 0; j < it->first.size(); j++) {
      buffer[bufpt++] = it->first[j];
    }
    buffer[bufpt++] = '\0';

    char flt[sizeof(float)];
    float val = it->second;
    memcpy(&flt, &val, sizeof(float));
    for (unsigned int j = 0; j < sizeof(float); j++) {
      buffer[bufpt++] = flt[j];
    }
  }

  * size = totalSize;
  return buffer;
}

void deserialize(vector<pair<string, float>> &restore,
                 char * buffer, unsigned int size) {
  bool term = false;
  string curStr;
  char flt[sizeof(float)];
  unsigned int fltptr = 0;

  for (unsigned int j = 0; j < size; j++) {
    const char c = buffer[j];
    if (!term) {  // Read in string before first appearence of '\0'
      if (c != '\0') {
        curStr.push_back(c);
      } else {
        term = true;
      }
    } else {  // Read in float
      if (fltptr < sizeof(float)) {
        flt[fltptr] = buffer[j];
        fltptr++;
      } else {
        float number;
        memcpy(&number, &flt, sizeof(float));
        restore.push_back(make_pair(curStr, number));

        // Reset vars
        term = false;
        fltptr = 0;
        curStr.clear();
        curStr.push_back(c);
      }
    }
  }
}
