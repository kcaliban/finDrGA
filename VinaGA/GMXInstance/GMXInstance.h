#ifndef GMXINST
#define GMXINST
#include <string>
#include <cstdlib>
#include <iostream>
class GMXInstance {
  // Problem: forcefield folder needs to be in working dir and mdppath
  public:
    GMXInstance(const char * ligand1, const char * gromacsPath1,
                const char * workDir1, const char * forcefield1,
                const char * water1, const char * bt1,
                float boxsize1, const char * mdpPath1) {
      ligand = ligand1;
      gromacsPath = gromacsPath1;
      workDir = workDir1;
      forcefield = forcefield1;
      water = water1;
      bt = bt1;
      boxsize = boxsize1;
      mdpPath = mdpPath1; // Where all parameter files are
    }

    // Prepares a PDB for MD simulation
    void preparePDB();
    // Run MD simulation
    void runMD();
    // Returns path to clustered MD simulation (after clstering, that is)
    std::string clusteredMD();

  private:
    std::string ligand;
    std::string workDir;
    std::string forcefield;
    std::string water;
    std::string bt;
    std::string gromacsPath;
    std::string mdpPath;
    float boxsize;
};

#endif
