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
  return std::get<2>(internalMap[FASTASEQ]);
}

void PoolMGR::addElement(std::string FASTASEQ) {
  if (internalMap.count(FASTASEQ) == 0) {
    // Add object to map
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
  command.append(" 2>/dev/null");
  int success = system(command.c_str());
  if (success == -1) {
    std::cout << "Error creating directory for PDB file!\n"
              << "Sequence name: " << FASTASEQ << std::endl;
    exit(-1);
  }
  command.clear();
  // Run Pymol to create ligand
  command.append("pymol -kcQ -d \"fab ");
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
  success = system(command.c_str());
  if (success == -1) {
    std::cout << "Error creating PDB file from sequence!\n"
              << "Sequence name: " << FASTASEQ << std::endl;
    exit(-1);
  }
  // Add path to map
  std::get<0>(internalMap[FASTASEQ]) = workDir + "/" + FASTASEQ + "/" + FASTASEQ + ".pdb";
}

void PoolMGR::genMD(std::string FASTASEQ) {
  // TODO, for testing use the same as PDB
  std::get<1>(internalMap[FASTASEQ]) = std::get<0>(internalMap[FASTASEQ]);
}

void PoolMGR::genDock(std::string FASTASEQ) {
  VinaInstance vinaInstance(vinaPath.c_str(), pythonShPath.c_str(),
                            mgltoolstilitiesPath.c_str(),
                            workDir.c_str(),
                            receptor.c_str(),
                            std::get<1>(internalMap[FASTASEQ]).c_str());
  // std::cout << "Trying to generate config for: " << FASTASEQ << std::endl;
  vinaInstance.generateConf();
  // std::cout << "Generated config for: " << FASTASEQ << std::endl;
  vinaInstance.generatePDBQT();
  std::get<2>(internalMap[FASTASEQ]) =
    vinaInstance.calculateBindingAffinity(exhaustiveness, energy_range);
}