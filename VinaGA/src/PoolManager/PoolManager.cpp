#include "PoolManager.h"

void PoolMGR::printSeqAff() {
  for (auto it = internalMap.begin(); it != internalMap.end(); it++) {
    std::cout << it->first << ": " << std::get<2>(it->second) << std::endl;
  }
}

float PoolMGR::getAffinity(std::string FASTASEQ) {
  if (internalMap.count(FASTASEQ) == 0) {
    // Sequence is not in our pool
    genPDB(FASTASEQ);
    genMD(FASTASEQ);
    genDock(FASTASEQ);
  }
  float affinity;
  #pragma omp critical
  affinity = std::get<2>(internalMap.at(FASTASEQ));
  return affinity;
}

void PoolMGR::addElement(std::string FASTASEQ) {
  int count;
  #pragma omp critical
  count = internalMap.count(FASTASEQ);
  if (count == 0) {
    // Add object to map
    #pragma omp critical
    internalMap[FASTASEQ] = std::make_tuple("", "", 0, 0);
    // Generate PDB, MD and fitness function
    genPDB(FASTASEQ);
    genMD(FASTASEQ);
    genDock(FASTASEQ);
  }
}

void PoolMGR::genPDB(std::string FASTASEQ) {
  // Create directory
  std::string command = "mkdir ";
  command.append(workDir);
  command.append("/");
  command.append(FASTASEQ);
  command.append(" 2>/dev/null 1>&2");
  int success = system(command.c_str());
  if (success == -1) {
    std::cout << "Error creating directory for PDB file!\n"
              << "Sequence name: " << FASTASEQ << std::endl;
    exit(-1);
  }
  command.clear();
  // Run Pymol to create ligand
  command.append(pymolPath);
  command.append(" -kcQ -d \"fab ");
  command.append(FASTASEQ);
  command.append(", ");
  command.append(FASTASEQ);
  command.append(", ss=1;"); // Secondary structure
  command.append("save ");
  command.append(workDir);
  command.append("/");
  command.append(FASTASEQ);
  command.append("/");
  command.append(FASTASEQ);
  command.append(".pdb\"");
  command.append(" >/dev/null 2>&1");
  success = system(command.c_str());
  if (success == -1) {
    std::cout << "Error creating PDB file from sequence!\n"
              << "Sequence name: " << FASTASEQ << std::endl;
    exit(-1);
  }
  // Add path to map
  #pragma omp critical
  std::get<0>(internalMap[FASTASEQ]) = workDir + "/" + FASTASEQ + "/" + FASTASEQ + ".pdb";
}

void PoolMGR::genMD(std::string FASTASEQ) {
  GMXInstance gmxInstance(std::get<0>(internalMap[FASTASEQ]).c_str(),
                          gromacsPath.c_str(), pymolPath.c_str(),
                          (workDir + "/" + FASTASEQ).c_str(),
                          forcefield.c_str(), forcefieldPath.c_str(), water.c_str(),
                          boundingboxtype.c_str(), clustercutoff,
                          boxsize, mdpPath.c_str(), true, true);
  gmxInstance.preparePDB();
  gmxInstance.runMD();
  gmxInstance.clusteredMD();
  gmxInstance.extractTopCluster();

  #pragma omp critical
  std::get<1>(internalMap[FASTASEQ]) = workDir + "/" + FASTASEQ + "/topcluster.pdb";
}

void PoolMGR::genDock(std::string FASTASEQ) {
  VinaInstance vinaInstance(vinaPath.c_str(), pythonShPath.c_str(),
                            mgltoolstilitiesPath.c_str(),
                            pymolPath.c_str(),
                            workDir.c_str(),
                            receptor.c_str(),
                            std::get<1>(internalMap[FASTASEQ]).c_str(),
                            true, true);
  // std::cout << "Trying to generate config for: " << FASTASEQ << std::endl;
  vinaInstance.generateConf();
  // std::cout << "Generated config for: " << FASTASEQ << std::endl;
  vinaInstance.generatePDBQT();
  float affinity = vinaInstance.calculateBindingAffinity(exhaustiveness, energy_range);

  #pragma omp critical
  std::get<2>(internalMap[FASTASEQ]) = affinity;
}

void PoolMGR::update(std::vector<std::string> gen) {
  // Set number of times not used to minus one for all elements that have been used
  for (auto FASTASEQ : gen) {
    std::get<3>(internalMap.at(FASTASEQ)) = -1;
  }
  // Increase every number by one
  for (auto it = internalMap.begin(); it != internalMap.end(); it++) {
    std::get<3>(it->second) += 1;
  }
};

void PoolMGR::cleanUp(int x) {
  // Cleans all files of PDBs that have not been used for x generations
  for (auto it = internalMap.begin(); it != internalMap.end(); ) {
    if (std::get<3>(it->second) > x) {
      std::cout << "CLEANUP\t\tRemoving: " << it->first << std::endl;
      deleteElementData(it->first);
      internalMap.erase(it++);
    } else {
      it++;
    }
  }
}

void PoolMGR::deleteElementData(std::string FASTASEQ) {
  std::cout << "Has not been used for some generations: "
            << FASTASEQ << std::endl;
  // Simply remove the directory recursively
  std::string cmd;
  cmd.append("rm -rf ");
  cmd.append(workDir);
  cmd.append("/");
  cmd.append(FASTASEQ);
  int success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Error deleting directory in cleaning step of PoolMGR!\n"
              << "Directory: " << workDir << "/" << FASTASEQ << std::endl;
    exit(-1);
  }
}
