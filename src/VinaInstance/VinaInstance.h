/* Copyright 2019 Fabian Krause
 *
 * AutoDock Vina Interface class
 *
 * Provides functionality to prepare and execute an AutoDock Vina docking
 * using system() calls.
 *
*/
#ifndef SRC_VINAINSTANCE_VINAINSTANCE_H_
#define SRC_VINAINSTANCE_VINAINSTANCE_H_
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <regex>
#include <sstream>
#include <fstream>
#include <exception>
#include "../Info.h"
class VinaException : virtual public std::exception {
 public:
    std::string type;
    VinaException(const std::string msg1, const std::string file1) {
      msg = msg1;
      file = file1;
      std::string error;
      error.append("Error in VinaInstance!\n");
      error.append("File: ");
      error.append(file);
      error.append("\n");
      error.append("Message: ");
      error.append(msg);
      errorMsg = error;
      type = "";
    }

    VinaException(const std::string msg1,
                  const std::string file1,
                  const std::string type1) {
      msg = msg1;
      file = file1;
      std::string error;
      error.append("Error in VinaInstance!\n");
      error.append("File: ");
      error.append(file);
      error.append("\n");
      error.append("Message: ");
      error.append(msg);
      errorMsg = error;
      type = type1;
    }

    virtual const char * what() const throw() {
      return errorMsg.c_str();
    }

 private:
    std::string msg;
    std::string file;
    std::string errorMsg;
};

class VinaInstance {
 public:
    VinaInstance(const char * vinaPath1,
                 const char * receptor1,
                 const char * ligand1,
                 Info * info1) {
      vinaPath = vinaPath1;
      receptor = receptor1;
      ligand = ligand1;
      info = info1;
    }

    /* calculateBindingAffinity(exhaustiveness, energy_range):
     *
     * Does a docking with specified receptor and ligand using
     * exhaustiveness and energy_range settings.
     * Returns the best docking result.
     *
    */
    float calculateBindingAffinity(int, int);

 private:
    std::string vinaPath;
    std::string receptor;
    std::string ligand;
    Info * info;
};

#endif  // SRC_VINAINSTANCE_VINAINSTANCE_H_

