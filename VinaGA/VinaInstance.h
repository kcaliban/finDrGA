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
#include <filesystem>
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
    VinaInstance(const char *  vinaPath1, const char * pythonShPath1,
                  const char * mgltoolstilitiesPath1, const char * workDir1,
                  const char * receptor1, const char * ligand1) {
      vinaPath = vinaPath1;
      pythonShPath = pythonShPath1;
      mgltoolstilitiesPath = mgltoolstilitiesPath1;
      workDir = workDir1;
      receptor = receptor1;
      ligand = ligand1;

      char cmd[2000];
      strcpy(cmd, "mkdir ");
      strcat(cmd, workDir.c_str());
      strcat(cmd, " 2>/dev/null"); // Error is usually that path already exists
      system(cmd);
    }
    void generatePDBQT();

    void generateConf();

    float calculateBindingAffinity(int exhaustiveness, int energy_range);

  private:
    std::string vinaPath;
    std::string pythonShPath;
    std::string mgltoolstilitiesPath;
    std::string receptor;
    std::string ligand;
    std::string workDir;
};
#endif
