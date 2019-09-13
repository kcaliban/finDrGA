#ifndef GMXINST
#define GMXINST
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <exception>
class GMXException : public std::exception {
  public:
    std::string type;

    GMXException(const std::string msg1, const std::string file1) {
      msg = msg1;
      file = file1;
      type = "";
      std::string error;
      error.append("Error in GMXInstance!\n");
      error.append("File: ");
      error.append(file);
      error.append("\n");
      error.append("Message: ");
      error.append(msg);
      errorMsg = error;
    }

    GMXException(const std::string msg1, const std::string file1,
                 const std::string type1) {
      msg = msg1;
      file = file1;
      type = type1;
      std::string error;
      error.append("Error in GMXInstance!\n");
      error.append("File: ");
      error.append(file);
      error.append("\n");
      error.append("Message: ");
      error.append(msg);
      errorMsg = error;
    }

    const char * what () const throw () {
      return errorMsg.c_str();
    }

  private:
    std::string msg;
    std::string file;
    std::string errorMsg;
};

class GMXInstance {
  public:
    GMXInstance(const char * ligand1, const char * gromacsPath1,
                const char * pymolPath1,
                const char * workDir1, const char * forcefield1,
                const char * forcefieldPath1,
                const char * water1, const char * bt1,
                float clustercutoff1,
                float boxsize1, const char * mdpPath1,
                bool debug1, bool log1) {
      ligand = ligand1;
      gromacsPath = gromacsPath1;
      workDir = workDir1;
      forcefield = forcefield1;
      water = water1;
      bt = bt1;
      boxsize = boxsize1;
      clustercutoff = clustercutoff1;
      mdpPath = mdpPath1;
      forcefieldPath = forcefieldPath1;
      pymolPath = pymolPath1;
      debug = debug1;
      log = log1;
    }

    void preparePDB();
    void runMD();
    void clusteredMD();
    void extractTopCluster();
    std::string logStr();

  private:
    std::string ligand;
    std::string workDir;
    std::string forcefield;
    std::string pymolPath;
    std::string water;
    std::string bt;
    std::string gromacsPath;
    std::string mdpPath;
    std::string forcefieldPath;
    float boxsize;
    float clustercutoff;
    bool debug;
    bool log;

    void debugPrint(const char *);
    void errorPrint(const char *);
};

#endif
