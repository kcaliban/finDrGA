#ifndef GMXINST
#define GMXINST
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
class GMXInstance {
  // Problem: forcefield folder needs to be in working dir and mdppath
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
      mdpPath = mdpPath1; // Where all parameter files are
      forcefieldPath = forcefieldPath1;
      pymolPath = pymolPath1;
      debug = debug1;
      log = log1;
    }

    // Prepares a PDB for MD simulation
    void preparePDB();
    // Run MD simulation
    void runMD();
    // Returns path to clustered MD simulation (after clstering, that is)
    void clusteredMD();
    // Extract the top clusters
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
