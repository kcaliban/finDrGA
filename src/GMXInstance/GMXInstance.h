/* Copyright 2019 Fabian Krause
 *
 * GROMACS Interface class
 *
 * Provides functionality to prepare and run a molecular dynamics
 * simulation using system() calls.
 *
 * Preparation and settings follow the steps from the tutorial
 * "Lysozyme in Water" by Justin A. Lemkuhl, Ph.D.
 * http://www.mdtutorials.com/gmx/lysozyme/index.html
 *
*/
#ifndef SRC_GMXINSTANCE_GMXINSTANCE_H_
#define SRC_GMXINSTANCE_GMXINSTANCE_H_
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <exception>
#include "../Info.h"
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

    GMXException(const std::string msg1,
                 const std::string file1,
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

    const char * what() const throw() {
      return errorMsg.c_str();
    }

 private:
    std::string msg;
    std::string file;
    std::string errorMsg;
};

class GMXInstance {
 public:
    GMXInstance(const char * ligand1,
                const char * gromacsPath1,
                const char * pymolPath1,
                const char * workDir1,
                const char * forcefield1,
                const char * forcefieldPath1,
                const char * water1,
                const char * bt1,
                float clustercutoff1,
                float boxsize1,
                const char * mdpPath1,
                Info * info1) {
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
      info = info1;
    }

    /* preparePDB():
     * Prepares the ligand for molecular dynamics simulation by performing:
     * 1) Cleansing from crystal water
     * 2) Creation of a topology using specified force-field
     * 3) Solvating
     * 4) Ionizing
     * 5) Minimizing energy by a short MD specified in minim.mdp
     * 6) Equilibrating temperature and pressure doing short MDs specified
     *    in npt.mdp & nvt.mdp
     *
     * Relevant output: md_0_1.tpr
    */
    void preparePDB();
    /* runMD():
     * Runs a molecular dynamics simulation using the settings specified in
     * MD.mdp using the prepared PDB file from preparePDB()
     *
     * Relevant output: MD.pdb
    */
    void runMD();
    /* clusterMD():
     * Clusters the result of molecular dynamics simulation using gmx cluster
     *
     * Relevant output: clusters.pdb
    */
    void clusterMD();
    /* extractTopCluster():
     * Extracts the biggest cluster from clustering result
     *
     * Relevant output: topcluster.pdb
    */
    void extractTopCluster();
    void energyMinim();

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
    Info * info;

    /* logStr():
     * Returns command-line string to redirect stdout and stderr to log file
     */
    std::string logStr();
};

#endif  // SRC_GMXINSTANCE_GMXINSTANCE_H_
