/* Copyright 2019 Fabian Krause
 *
 * Info
 *
 * Manages logging and console output
*/
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

    /* infoMsg(str):
     *
     * Prints str to console in green and logfile according to settings
    */
    void infoMsg(std::string);
    /* errorMsg(str, fatal):
     *
     * Prints str to console in red and logfile according to settings,
     * exiting the program if fatal is true
    */
    void errorMsg(std::string, bool);

 private:
    bool log;
    bool console;
    std::ofstream * f;
};

#endif  // SRC_INFO_H_
