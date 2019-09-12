#ifndef POOLMGR
#define POOLMGR
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <vector>
#include "../VinaInstance/VinaInstance.h"
#include "../GMXInstance/GMXInstance.h"
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

    // Add PDB file, return value is FASTA seq
    std::string addElementPDB(std::string);
    // Adds an element if it does not exist already
    void addElement(std::string);
    // Read out affinity for given string
    float getAffinity(std::string);
    // Updates number of rounds not used; input: vector of FASTA seq.
    // in current Generation
    void update(std::vector<std::string>);
    // Deletes any element and their MD & PDB files that have not been used
    // for int generations
    void cleanUp(int);
    // Print all sequences and their affinity
    void printSeqAff();
    // The above as a string
    std::string toStr();
  private:
    // Exhaustiveness & Energy range
    int exhaustiveness;
    int energy_range;
    // Goal receptors
    std::vector<std::string> receptors;
    int nReceptors;
    // Working directory
    std::string workDir;
    // Directories required for AutoDock VINA
    std::string vinaPath;
    std::string pythonShPath;
    std::string mgltoolstilitiesPath;
    std::string pymolPath;
    // Directories required for GROMACS
    std::string gromacsPath;
    std::string mdpPath;
    std::string forcefield;
    std::string forcefieldPath;
    std::string water;
    std::string boundingboxtype;
    float boxsize;
    float clustercutoff;
    // Map of FASTA sequence and a tuple containing (if applicable)
    // path of PDB file, path of MD file, result of fitness function,
    // integer of number of rounds not accessed to
    std::unordered_map<std::string, std::tuple<std::string,
                                               std::string,
                                               std::vector<int>,
                                               int> > internalMap;

    // Generate PDB
    void genPDB(std::string);
    // Generate MD
    void genMD(std::string);
    // Get docking result
    void genDock(std::string);
    // Delete an element
    void deleteElementData(std::string);
    // Generate FASTA from file
    std::string PDBtoFASTA(std::string);
};

#endif
