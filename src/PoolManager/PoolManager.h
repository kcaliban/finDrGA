/* Copyright 2019 Fabian Krause
 *
 * PoolManager
 *
 * Manages the gene pool, i.e. all FASTA sequences, their PDB files
 * and their MD as well as docking results.
 *
 * Can automatically delete unused files after certain number of generations
 * of non-usage.
*/
#ifndef SRC_POOLMANAGER_POOLMANAGER_H_
#define SRC_POOLMANAGER_POOLMANAGER_H_
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <exception>
#include <string>
#include <utility>
#include "../VinaInstance/VinaInstance.h"
#include "../GMXInstance/GMXInstance.h"
#include "../Serialization/Serialization.h"
#include "../Communication.h"
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

    virtual const char * what() const throw() {
      return errorMsg.c_str();
    }

 private:
    std::string msg;
    std::string file;
    std::string errorMsg;
};

class PoolMGR {
 public:
    PoolMGR(const char * workDir1,
            const char * vinaPath1,
            const char * pythonShPath1,
            const char * mgltoolstilitiesPath1,
            const char * pymolPath1,
            std::vector<std::string> receptors1,
            int exhaustiveness1,
            int energy_range1,
            const char * gromacsPath1,
            const char * mdpPath1,
            const char * forcefield1,
            const char * forcefieldPath1,
            const char * water1,
            const char * boundingboxtype1,
            float boxsize1,
            float clustercutoff1,
            Info * info1,
            bool pymolgen1) {
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
      info = info1;
      pymolgen = pymolgen1;
    }

    /* addElementPDB(path):
     *
     * Adds an existing .pdb file to the manager by generating its FASTA
     * sequence, coyping it into the working directory and performing
     * a Docking and MD
     *
     * Returns FASTASEQ if MD and Docking successful, else empty string
    */
    std::string addElementPDB(std::string, bool);
    /* addElemennt(FASTA):
     *
     * Creates a PDB in alpha-helical structure using pymol from FASTA sequence
     * and adds it to the manager
    */
    void addElement(std::string);
    /* getAffinity(FASTA):
     *
     * Returns the calculated affinity
    */
    float getAffinity(std::string);
    /* update(vector of FASTAs):
     *
     * Updates number of rounds unused for internal gene pool
    */
    void update(std::vector<std::string>);
    /* cleanUp(n):
     *
     * Removes PDBs not used for n rounds
    */
    void cleanUp(int);
    /* preparePDBQT(FASTA):
     *
     * Prepares a ligand for docking for each receptor
    */
    void preparePDBQT(std::string);
    /* toStr():
     *
     * Returns a string containing every individual in the gene pool and their
     * docking results to each receptor
    */
    std::string toStr();
    /* PDBtoFASTA(path):
     *
     * Returns FASTA sequence of specified file
    */
    std::string PDBtoFASTA(std::string);
    /* addElementsFromFASTAs(FASTAs, world_size):
     *
     * Adds elements to the pool from FASTA sequences, distributing amongst
     * nodes according to available threads
    */
    std::vector<std::string> addElementsFromFASTAs(std::vector<std::string>&,
                                                   int);
    /* addElementsFromPDBs(PDB paths, world_size):
     *
     * Adds elements to the pool from PDB file paths, distributing amongst
     * nodes according to available threads
    */
    std::vector<std::string> addElementsFromPDBs(std::vector<std::string>&,
                                                 int);
    /* getFASTAS(PDB paths):
     *
     * Collect FASTA sequences for given PDB file path vector
    */
    std::vector<std::string> getFASTAS(std::vector<std::string> &);

 private:
    int exhaustiveness;
    int energy_range;
    std::vector<std::string> receptors;
    std::vector<unsigned int> jobDistribution;
    int nReceptors;
    std::string workDir;
    std::string vinaPath;
    std::string pythonShPath;
    std::string mgltoolstilitiesPath;
    std::string pymolPath;
    Info * info;

    std::string gromacsPath;
    std::string mdpPath;
    std::string forcefield;
    std::string forcefieldPath;
    std::string water;
    std::string boundingboxtype;
    float boxsize;
    float clustercutoff;
    std::unordered_map<std::string,
                       std::tuple<std::string,
                                  std::string,
                                  float,
                                  int> > internalMap;
    bool pymolgen;

    /* genPDB(FASTA):
     *
     * Generates a PDB in alpha helical structure using pymol fab
    */
    void genPDB(std::string);
    /* genMD(FASTA):
     *
     * Performs a molecular dynamics simulation for
     * specified FASTA sequence (if already in manager)
     * using GMXInstance
    */
    void genMD(std::string);
    /* genDock(FASTA):
     *
     * Performs a docking for
     * specified FASTA sequence (if already in manager)
     * using VinaInstance
    */
    void genDock(std::string);
    /* deleteElementData(FASTA)
     *
     * Deletes all files for specified FASTA sequence
     *
    */
    void deleteElementData(std::string);
    /* addElementsFromFiles(File paths, world_size):
     *
     * Used by addElementsFromPDBs and addElementsFromFASTAs to distribute
     * docking and MD to computing nodes and collect the results
     *
    */
    std::vector<std::string> addElementsFromFiles(std::vector<std::string>&,
                                                  int);
};

#endif  //  SRC_POOLMANAGER_POOLMANAGER_H_
