#ifndef VININS
#define VININS
/* Interface to run Vina and read simulation results */
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
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

    VinaException(const std::string msg1, const std::string file1, std::string type1) {
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

    virtual const char * what () const throw () {
      return errorMsg.c_str();
    }

  private:
    std::string msg;
    std::string file;
    std::string errorMsg;
};

class VinaInstance {
  public:
    /* Docking instance using AutoDock Vina
     * Paths are self-explanatory, but beware:
     *    ALWAYS WITHOUT '/' AT THE END
     *    NEVER STARTING WITH '~'
     *
     * workDir ~ working directory, programm will create subdirectory
     *           with name of ligand
     * receptor, ligand ~ absolute path to files
     */
    VinaInstance(const char *  vinaPath1, const char * receptor1,
                 const char * ligand1, Info * info1) {
      vinaPath = vinaPath1;
      receptor = receptor1;
      ligand = ligand1;
      info = info1;
    }

    float calculateBindingAffinity(int exhaustiveness, int energy_range);

  private:
    std::string vinaPath;
    std::string receptor;
    std::string ligand;
    Info * info;
};

#endif
