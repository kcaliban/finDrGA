#ifndef GMXINST
#define GMXINST
#include <string>
#include <cstdlib>
#include <iostream>
class GMXInstance {
  public:
    GMXInstance(const char * ligand1, const char * gromacsPath1,
                const char * workDir1, const char * forcefield1,
                const char * water1, const char * bt1) {
      ligand = ligand1;
      gromacsPath = gromacsPath1;
      workDir = workDir1;
      forcefield = forcefield1;
      water = water1;
      bt = bt1;
    }

    // Prepares a PDB for MD simulation
    void preparePDB();
    // Returns path to clustered MD simulation
    std::string MD();
  private:
    std::string ligand;
    std::string workDir;
    std::string forcefield;
    std::string water;
    std::string bt;
    std::string gromacsPath;
};

#endif
