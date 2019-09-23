/* Copyright 2019 Fabian Krause */
#include "PoolManager.h"

void PoolMGR::genMD(std::string) {
  return;
}

void PoolMGR::preparePDBQT(std::string ligand) {
  info->infoMsg("(POOLMGR) Preparing PDBQT of ligand: " + ligand);
  // Generate a PDBQT
  std::string command;
  command.append(pythonShPath);
  command.append(" ");
  command.append(mgltoolstilitiesPath);
  command.append("/prepare_ligand4.py -l ");
  command.append(ligand);
  command.append(" -Z -A bonds_hydrogens -U nphs -o ");
  command.append(ligand);
  command.append("qt >/dev/null");
  int success = system(command.c_str());
  if (success != 0) {
    throw VinaException("Could not generate pdbqt file for ligand",
                        ligand,
                        "PQT");
  }
}

std::string PoolMGR::PDBtoFASTA(std::string filename) {
  info->infoMsg("(POOLMGR) Getting FASTA from PDB: " + filename);
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw PoolManagerException(
          "Could not open PDB file to convert to FASTA sequence",
          filename);
  }
  std::unordered_map<std::string, std::string> AA({
    {"ALA", "A"}, {"ARG", "R"}, {"ASN", "N"}, {"ASP", "D"}, {"CYS", "C"},
    {"GLU", "E"}, {"GLN", "Q"}, {"GLY", "G"}, {"HIS", "H"}, {"ILE", "I"},
    {"LEU", "L"}, {"LYS", "K"}, {"MET", "M"}, {"PHE", "F"}, {"PRO", "P"},
    {"SER", "S"}, {"THR", "T"}, {"TRP", "W"}, {"TYR", "Y"}, {"VAL", "V"}
  });

  std::string FASTA;
  std::string line;
  int previd = -1;
  while (getline(file, line)) {
    if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
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

std::vector<std::string> PoolMGR::addElementsFromPDBs(std::vector<std::string> &files,
                                  int world_size) {
  // Prepare internal map and directory structure
  std::vector<std::string> newFiles;
  for (auto file : files) {
    // Get FASTA sequence
    std::string FASTASEQ = PDBtoFASTA(file);
    // If error in FASTA generation or already in pool, continue
    if (FASTASEQ.empty() || internalMap.count(FASTASEQ) != 0) {continue;}
    internalMap[FASTASEQ] = std::make_tuple(workDir + "/" + FASTASEQ + "/" +
                                            FASTASEQ + ".pdb",
                                            workDir + "/" + FASTASEQ + "/" +
                                            FASTASEQ + ".pdb",
                                            10.0f, 0);
    // Make required directory
    std::string command = "mkdir -p ";
    command.append(workDir);
    command.append("/");
    command.append(FASTASEQ);
    command.append(" 2>/dev/null 1>&2");
    int success = system(command.c_str());
    if (success != 0) {
      PoolManagerException("Could not create directory for PDB file", FASTASEQ);
    }
    command.clear();
    // Copy file
    command.append("cp -f ");
    command.append(file);
    command.append(" ");
    command.append(workDir);
    command.append("/");
    command.append(FASTASEQ);
    command.append("/");
    command.append(FASTASEQ);
    command.append(".pdb");
    success = system(command.c_str());
    if (success != 0) {
      throw PoolManagerException("Could not copy PDB file", file);
    }
    command.clear();
    newFiles.push_back(workDir + "/" + FASTASEQ + "/" + FASTASEQ + ".pdb");
  }
  // Do MD and get docking results from PoolWorker
  // Spread files to available subprocesses
  // Exclude main process from calculation hence world_size - 1
  unsigned int perNode = ceil((float) newFiles.size() / (float) (world_size - 1));
  std::vector<std::vector<std::string>> buckets;
  for (int i = 0; i < world_size - 1; i++) {
    std::vector<std::string> bucket;
    for (unsigned int j = 0; j < perNode; j++) {
      if (i * perNode + j < newFiles.size()) {
        bucket.push_back(newFiles.at(i * perNode + j));
      } else {
        break;
      }
    }
    buckets.push_back(bucket);
  }
  // Send buckets to subprocesses
  for (int i = 1; i < world_size; i++) {
    info->infoMsg("Bossman is sending his work...");
    unsigned int bucketSize;
    char * tmp = serialize(buckets.at(i - 1), &bucketSize);
    MPI_Send(&bucketSize, 1, MPI_INT, i, SENDFILESSIZE, MPI_COMM_WORLD);
    MPI_Send(&tmp[0], bucketSize, MPI_BYTE, i, SENDFILESCONT, MPI_COMM_WORLD);
    free(tmp);
  }
  // Get results of each bucket
  std::vector<std::vector<std::pair<std::string, float>>> bucketResults;
  for (int i = 1; i < world_size; i++) {
    std::vector<std::pair<std::string, float>> results;
    unsigned int resultsSize;

    MPI_Recv(&resultsSize, 1, MPI_INT, i, SENDAFFINSIZE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    char * tmp = new char[resultsSize];
    MPI_Recv(&tmp[0], resultsSize, MPI_BYTE, i, SENDAFFINCONT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    deserialize(results, tmp, resultsSize);
    bucketResults.push_back(results);
    free(tmp);
  }
  // Add the results to map and return FASTA sequences of add results
  std::vector<std::string> returnVal;
  for (int i = 1; i < world_size; i++) {
    for (unsigned int j = 0; j < bucketResults.at(i - 1).size(); j++) {
      std::string path = bucketResults.at(i - 1).at(j).first;
      size_t lastSlash = path.find_last_of("/");
      size_t dot = path.find_last_of(".");
      std::string fasta = path.substr(lastSlash + 1,
                                      path.size() - lastSlash
                                      - (path.size() - dot + 1));
      returnVal.push_back(fasta);
      float aff = bucketResults.at(i - 1).at(j).second;
      // Next MD continues where previous MD has left off
      std::get<1>(internalMap[fasta]) = workDir + "/"
                                        + fasta + "/topcluster.pdb";
      std::get<2>(internalMap[fasta]) = aff;
    }
  }
  return returnVal;
}

std::vector<std::string> PoolMGR::addElementsFromFASTAs(
                                    std::vector<std::string> &fastas,
                                    int world_size) {
  // Prepare files
  for (auto i : fastas) {
    if (internalMap.count(i) == 0) {
      internalMap[i] = std::make_tuple("", "", 10.0f, 0);
      genPDB(i);
    }
  }
  std::vector<std::string> newFiles;
  for (auto i : fastas) {
    newFiles.push_back(std::get<1>(internalMap[i]));
  }
  // Do MD and get docking results from PoolWorker
  // Spread files to available subprocesses
  // Exclude main process from calculation hence world_size - 1
  unsigned int perNode = ceil((float) newFiles.size() / (float) (world_size - 1));
  std::vector<std::vector<std::string>> buckets;
  for (int i = 0; i < world_size - 1; i++) {
    std::vector<std::string> bucket;
    for (unsigned int j = 0; j < perNode; j++) {
      if (i * perNode + j < newFiles.size()) {
        bucket.push_back(newFiles.at(i * perNode + j));
      } else {
        break;
      }
    }
    buckets.push_back(bucket);
  }
  // Send buckets to subprocesses
  for (int i = 1; i < world_size; i++) {
    info->infoMsg("Bossman is sending his work...");
    unsigned int bucketSize;
    char * tmp = serialize(buckets.at(i - 1), &bucketSize);
    MPI_Send(&bucketSize, 1, MPI_INT, i, SENDFILESSIZE, MPI_COMM_WORLD);
    MPI_Send(&tmp[0], bucketSize, MPI_BYTE, i, SENDFILESCONT, MPI_COMM_WORLD);
    free(tmp);
  }
  // Get results of each bucket
  std::vector<std::vector<std::pair<std::string, float>>> bucketResults;
  for (int i = 1; i < world_size; i++) {
    std::vector<std::pair<std::string, float>> results;
    unsigned int resultsSize;

    MPI_Recv(&resultsSize, 1, MPI_INT, i, SENDAFFINSIZE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    char * tmp = new char[resultsSize];
    MPI_Recv(&tmp[0], resultsSize, MPI_BYTE, i, SENDAFFINCONT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    deserialize(results, tmp, resultsSize);
    bucketResults.push_back(results);
    free(tmp);
  }
  // Add the results to map and return FASTA sequences of added results
  std::vector<std::string> returnVal;
  for (int i = 1; i < world_size; i++) {
    for (unsigned int j = 0; j < bucketResults.at(i - 1).size(); j++) {
      std::string path = bucketResults.at(i - 1).at(j).first;
      size_t lastSlash = path.find_last_of("/");
      size_t dot = path.find_last_of(".");
      std::string fasta = path.substr(lastSlash + 1,
                                      path.size() - lastSlash
                                      - (path.size() - dot + 1));
      returnVal.push_back(fasta);
      float aff = bucketResults.at(i - 1).at(j).second;
      std::get<2>(internalMap[fasta]) = aff;
    }
  }
  return returnVal;
}

std::string PoolMGR::toStr() {
  std::string returnStr;
  returnStr.append("[");
  for (auto it = internalMap.begin(); it != internalMap.end(); it++) {
    returnStr.append(it->first);
    returnStr.append(": ");
    returnStr.append(std::to_string(std::get<2>(it->second)));
    returnStr.append("; ");
  }
  returnStr.pop_back(); returnStr.pop_back();
  returnStr.append("]");

  return returnStr;
}

float PoolMGR::getAffinity(std::string FASTASEQ) {
  return std::get<2>(internalMap.at(FASTASEQ));
}

void PoolMGR::genPDB(std::string FASTASEQ) {
  // Create directory
  std::string command = "mkdir -p ";
  command.append(workDir);
  command.append("/");
  command.append(FASTASEQ);
  command.append(" 2>/dev/null 1>&2");
  int success = system(command.c_str());
  if (success != 0) {
    throw PoolManagerException("Could not create directory for PDB file",
                               FASTASEQ);
  }
  command.clear();
  // Run Pymol to create ligand
  command.append(pymolPath);
  command.append(" -kcQ -d \"fab ");
  command.append(FASTASEQ);
  command.append(", ");
  command.append(FASTASEQ);
  command.append(", ss=1;");  // Secondary structure
  command.append("save ");
  command.append(workDir);
  command.append("/");
  command.append(FASTASEQ);
  command.append("/");
  command.append(FASTASEQ);
  command.append(".pdb\"");
  command.append(" >/dev/null 2>&1");
  success = system(command.c_str());
  if (success != 0) {
    throw PoolManagerException("Could not create PDB file", FASTASEQ);
  }
  // Add path to map
  #pragma omp critical
  std::get<0>(internalMap[FASTASEQ]) =
                          workDir + "/" + FASTASEQ + "/" + FASTASEQ + ".pdb";
  // This will later be overwritten with MD'ed file
  std::get<1>(internalMap[FASTASEQ]) =
                          workDir + "/" + FASTASEQ + "/" + FASTASEQ + ".pdb";
}


void PoolMGR::update(std::vector<std::string> gen) {
  for (auto FASTASEQ : gen) {
    std::get<3>(internalMap.at(FASTASEQ)) = -1;
  }
  // Increase every number by one
  for (auto it = internalMap.begin(); it != internalMap.end(); it++) {
    std::get<3>(it->second) += 1;
  }
}

void PoolMGR::cleanUp(int x) {
  // Cleans all files of PDBs that have not been used for x generations
  for (auto it = internalMap.begin(); it != internalMap.end(); ) {
    if (std::get<3>(it->second) > x) {
      deleteElementData(it->first);
      free(&std::get<2>(it->second));  // Free up heap space
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
  std::string command;
  command.append("rm -rf ");
  command.append(workDir);
  command.append("/");
  command.append(FASTASEQ);
  int success = system(command.c_str());
  if (success != 0) {
    throw PoolManagerException("Could not clean unused PDB files",
              workDir + "/" + FASTASEQ);
  }
}
