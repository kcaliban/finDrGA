/* Copyright 2019 Fabian Krause */
#include "Info.h"

void Info::infoMsg(std::string msg) {
  if (console) {
    std::cout << "\033[1;32mINFO: " + msg << "\033[0m" << std::endl;
  }
  if (log) {
    *f << "INFO: " + msg << std::endl;
  }
}

void Info::errorMsg(std::string msg, bool fatal = false) {
  std::cout << "\033[1;31mERROR: " + msg << "\033[0m" << std::endl;
  if (log) {
    *f << "ERROR: " + msg << std::endl;
  }
  if (fatal) {
    exit(-1);
  }
}
