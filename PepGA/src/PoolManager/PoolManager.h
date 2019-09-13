#ifndef POOLMGR
#define POOLMGR
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <vector>
#include "../VinaInstance/VinaInstance.h"
#include "../GMXInstance/GMXInstance.h"
#include <exception>
class PoolManagerException : virtual public std::exception {
  public:
    PoolManagerException(const std::string msg1, const std::string file1) {
      msg = msg1;
      file = file1;
      std::string error;
      error.append("Error in PoolManager!\n");
      error.append("File: ");
      error.append(file);
      error.append("\n");
      error.append("Message: ");
      error.append(msg);
      errorMsg = error;
    }

    virtual const char * what () const throw () {
      return errorMsg.c_str();
    }

  private:
    std::string msg;
    std::string file;
    std::string errorMsg;
};

class PoolMGR {
  /* Pool Manager: Keeps track of all PDB files and their corresponding
   *               FASTA sequence.
   *               Also keeps track of their MD simulation and docking
   *               results, if applicable.
   *
   *               If a certain FASTA sequence has not been used for
   *               X
   *               generations, delete the files and entry in PoolMGR
   */
  public:
    PoolMGR(const char * workDir1, const char * vinaPath1,
            const char * pythonShPath1, const char * mgltoolstilitiesPath1,
            const char * pymolPath1,
            std::vector<std::string> receptors1,
            int exhaustiveness1,
            int energy_range1,
            const char * gromacsPath1, const char * mdpPath1,
            const char * forcefield1, const char * forcefieldPath1,
            const char * water1, const char * boundingboxtype1,
            float boxsize1, float clustercutoff1) {
      workDir = workDir1;
      receptors = receptors1;
      nReceptors = receptors1.size();
      vinaPath = vinaPath1;
      pythonShPath = pythonShPath1;
      pymolPath = pymolPath1;
      mgltoolstilitiesPath = mgltoolstilitiesPath1;
      exhaustiveness = exhaustiveness1;
      energy_range = energy_range1;
      gromacsPath = gromacsPath1;
      mdpPath = mdpPath1;
      forcefield = forcefield1;
      forcefieldPath = forcefieldPath1;
      water = water1;
      boundingboxtype = boundingboxtype1;
      boxsize = boxsize1;
      clustercutoff = clustercutoff1;
    };

    std::string addElementPDB(std::string);
    void addElement(std::string);
    float getAffinity(std::string);
    void update(std::vector<std::string>);
    void cleanUp(int);
    void printSeqAff();
    std::string toStr();
  private:
    int exhaustiveness;
    int energy_range;
    std::vector<std::string> receptors;
    int nReceptors;
    std::string workDir;
    std::string vinaPath;
    std::string pythonShPath;
    std::string mgltoolstilitiesPath;
    std::string pymolPath;

    std::string gromacsPath;
    std::string mdpPath;
    std::string forcefield;
    std::string forcefieldPath;
    std::string water;
    std::string boundingboxtype;
    float boxsize;
    float clustercutoff;
    std::unordered_map<std::string, std::tuple<std::string,
                                               std::string,
                                               std::vector<std::pair<std::string, float>>,
                                               int> > internalMap;

    void genPDB(std::string);
    void genMD(std::string);
    void genDock(std::string);
    void deleteElementData(std::string);
    std::string PDBtoFASTA(std::string);
};

#endif
