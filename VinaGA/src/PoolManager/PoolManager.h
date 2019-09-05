#ifndef POOLMGR
#define POOLMGR
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <vector>
#include "../VinaInstance/VinaInstance.h"
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
            const char * receptor1,
            int exhaustiveness1,
            int energy_range1) {
      workDir = workDir1;
      receptor = receptor1;
      vinaPath = vinaPath1;
      pythonShPath = pythonShPath1;
      mgltoolstilitiesPath = mgltoolstilitiesPath1;
      exhaustiveness = exhaustiveness1;
      energy_range = energy_range1;
    };

    // Add PDB file
    void addElementPDB(std::string);
    // Adds an element if it does not exist already
    void addElement(std::string);
    // Read out affinity for given string
    float getAffinity(std::string);
    // Updates number of rounds not used; input: vector of FASTA seq.
    // in current Generation
    void update(std::vector<std::string>);
    // Deletes any element and their MD & PDB files that have not been used
    // for 3 consecutive generations
    void cleanUp(int);
    // Print all sequences and their affinity
    void printSeqAff();
  private:
    // Exhaustiveness & Energy range
    int exhaustiveness;
    int energy_range;
    // Goal receptor clustered MD pdb
    std::string receptor;
    // Working directory
    std::string workDir;
    // Directories required for AutoDock VINA
    std::string vinaPath;
    std::string pythonShPath;
    std::string mgltoolstilitiesPath;
    // Map of FASTA sequence and a tuple containing (if applicable)
    // path of PDB file, path of MD file, result of fitness function,
    // integer of number of rounds not accessed to
    std::unordered_map<std::string, std::tuple<std::string,
                                               std::string,
                                               float,
                                               int> > internalMap;

    // Generate PDB
    void genPDB(std::string);
    // Generate MD
    void genMD(std::string);
    // Get docking result
    void genDock(std::string);
};

#endif
