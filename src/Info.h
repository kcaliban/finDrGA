/* Copyright 2019 Fabian Krause */
#ifndef SRC_INFO_H_
#define SRC_INFO_H_
#include <string>
#include <iostream>
#include <fstream>
/* Logging, information and errors */
class Info {
 public:
    Info(bool log1, bool console1, std::string logfile1) {
      log = log1;
      console = console1;

      if (log) {
        f = new std::ofstream(logfile1, std::ios::out | std::ios::app);
      }
    }

    ~Info() {
      f->close();
      free(f);
    }

    void infoMsg(std::string);
    void errorMsg(std::string, bool);
 private:
    bool log;
    bool console;
    std::ofstream * f;
};
#endif  // SRC_INFO_H_
