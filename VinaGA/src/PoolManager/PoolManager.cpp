#include "PoolManager.h"

void PoolMGR::printSeqAff() {
  for (auto it = internalMap.begin(); it != internalMap.end(); it++) {
    std::cout << it->first << ": " << std::get<2>(it->second) << std::endl;
  }
}

std::string PoolMGR::PDBtoFASTA(std::string filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cout << "Error: Could not open PDB file to convert to FASTA sequence"
              << std::endl;
    std::cout << "Filename: " << filename << std::endl;
    exit(-1);
  }
  std::unordered_map<std::string, std::string> AA({
    {"ALA","A"},{"ARG","R"},{"ASN","N"},{"ASP","D"},{"CYS","C"},{"GLU","E"},{"GLN","Q"},{"GLY","G"},{"HIS","H"},
     {"ILE","I"},{"LEU","L"},{"LYS","K"},{"MET","M"},{"PHE","F"},{"PRO","P"},{"SER","S"},{"THR","T"},{"TRP","W"},
     {"TYR","Y"},{"VAL","V"}
  });

  std::string FASTA;
  std::string line;
  int previd = -1;
  while (getline(file, line)) {
    if (line.substr(0, 4) == "ATOM" or line.substr(0, 6) == "HETATM") {
      std::string aa = line.substr(17, 3);
      int id = stoi(line.substr(22, 4));
      if (id != previd) {
        FASTA.append(AA[aa]);
      }
      previd = id;
    }
  }
  return FASTA;
}

void PoolMGR::addElementPDB(std::string file) {
  // Get FASTA sequence
  std::string FASTASEQ = PDBtoFASTA(file);
  // If FASTA is already in pool we can return
  int count;
  #pragma omp critical
  count = internalMap.count(FASTASEQ);
  if (count != 0) { return;}
  // Make required directory
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
  // Move file
  command.append("mv ");
  command.append(file);
  command.append(" ");
  command.append(workDir);
  command.append("/");
  command.append(FASTASEQ);
  command.append("/");
  command.append(FASTASEQ);
  command.append(".pdb");
  success = system(command.c_str());
  if (success == -1) {
    std::cout << "Error moving PDB file!"
              << "File: " << file << std::endl;
    exit(-1);
  }
  command.clear();
  // Add file to internal map, do MD, dock
  #pragma omp critical
  internalMap[FASTASEQ] = std::make_tuple(FASTASEQ,
                                          workDir + "/" + FASTASEQ + "/" +
                                          FASTASEQ + ".pdb", 0, 0);
  genMD(FASTASEQ);
  genDock(FASTASEQ);
}

std::string PoolMGR::toStr() {
  std::string returnStr;
  returnStr.append("[");
  for (auto it = internalMap.begin(); it != internalMap.end(); it++) {
    returnStr.append(it->first);
    returnStr.append(": ");
    returnStr.append(std::to_string(std::get<2>(it->second)));
    returnStr.append(", ");
  }
  returnStr.pop_back(); returnStr.pop_back();
  returnStr.append("]");

  return returnStr;
}

float PoolMGR::getAffinity(std::string FASTASEQ) {
  /* Obsolete, we always add the element before we access its affinity
  if (internalMap.count(FASTASEQ) == 0) {
    std::cout << "Affinity requested even tho not in pool!" << std::endl;
    // Sequence is not in our pool
    genPDB(FASTASEQ);
    genMD(FASTASEQ);
    genDock(FASTASEQ);
  }
  */
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
  } else {
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
  std::string prevMD;
  #pragma omp critical
  prevMD = std::get<1>(internalMap[FASTASEQ]);
  std::string mdinput;
  // If there has not been an MD yet, we simply take the original pdb
  if (prevMD == "") {
    #pragma omp critical
    mdinput = std::get<0>(internalMap[FASTASEQ]);
  }
  // If there has, use it to do a further MD
  else {
    #pragma omp critical
    mdinput = std::get<1>(internalMap[FASTASEQ]);
  }
  GMXInstance gmxInstance(mdinput.c_str(),
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
