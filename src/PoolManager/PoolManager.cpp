/* Copyright 2019 Fabian Krause */
#include "PoolManager.h"

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

std::vector<unsigned int> distributeJobs(int * threads, int world_size) {
  // threads ~ number of avail threads per worker_i, i in (1, world_size - 1)
  // Simple approach: every computer gets as many as threads as possible
  std::vector<unsigned int> jobDistr;
  for (int j = 0; j < world_size - 1; j++) {
    jobDistr.push_back(threads[j]);
  }
  return jobDistr;
}

std::vector<std::string> PoolMGR::getFASTAS(std::vector<std::string> &files) {
  std::vector<std::string> FASTAS;
  for (auto file : files) {
    FASTAS.push_back(PDBtoFASTA(file));
  }
  return FASTAS;
}

std::vector<std::string> PoolMGR::addElementsFromFiles(std::vector<std::string>
                                                       &files, int world_size){
  info->infoMsg("Total number of affinities to be calculated: "
                 + std::to_string(files.size()));
  // Do MD and get docking results from PoolWorker
  int threadsPerWorker[world_size - 1];
  // Get the number of available threads for each worker to distribute
  // accordingly
  unsigned int allThreads = 0;
  for (int i = 1; i < world_size; i++) {
    unsigned int availThreads = 0;
    MPI_Recv(&availThreads, 1, MPI_INT, i,
             SENDNMTHREADS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    threadsPerWorker[i - 1] = availThreads;
    allThreads += availThreads;
  }
  info->infoMsg("Total number of available threads: "
                + std::to_string(allThreads));
  jobDistribution = distributeJobs(threadsPerWorker, world_size);
  std::vector<std::vector<std::string>> buckets;
  unsigned int vecPos = 0;
  for (int i = 0; i < world_size - 1; i++) {
    std::vector<std::string> bucket;

    for (unsigned int j = 0; j < jobDistribution.at(i); j++) {
      if (vecPos < files.size()) {
        bucket.push_back(files.at(vecPos++));
      } else {
        break;
      }
    }
    buckets.push_back(bucket);
  }
  // Send buckets to subprocesses
  info->infoMsg("Master is sending his work...");
  MPI_Request requests[world_size - 1];
  unsigned int bucketSizes[world_size - 1];
  char * bucketBin[world_size - 1];
  // Send size in bytes
  for (int i = 1; i < world_size; i++) {
    bucketBin[i - 1] = serialize(buckets.at(i - 1), &bucketSizes[i - 1]);
    MPI_Isend(&bucketSizes[i - 1], 1, MPI_INT, i,
              SENDFILESSIZE, MPI_COMM_WORLD, &requests[i - 1]);
  }
  MPI_Waitall(world_size - 1, requests, MPI_STATUS_IGNORE);
  // Send binary data
  for (int i = 1; i < world_size; i++) {
    MPI_Isend(&bucketBin[i - 1][0], bucketSizes[i - 1], MPI_BYTE, i,
              SENDFILESCONT, MPI_COMM_WORLD, &requests[i - 1]);
  }
  MPI_Waitall(world_size - 1, requests, MPI_STATUS_IGNORE);
  // Free memory
  for (int i = 1; i < world_size; i++) {
    free(bucketBin[i - 1]);
  }
  // Get results of each bucket
  std::vector<std::vector<std::pair<std::string, float>>> bucketResults;
  info->infoMsg("Master waiting for all results...");
  // Get filesize from each
  unsigned int resSize[world_size - 1];
  for (int i = 1; i < world_size; i++) {
    MPI_Irecv(&resSize[i - 1], 1, MPI_INT, i,
              SENDAFFINSIZE, MPI_COMM_WORLD, &requests[i - 1]);
  }
  MPI_Waitall(world_size - 1, requests, MPI_STATUS_IGNORE);
  info->infoMsg("Master got all filesizes");
  // Get contents from each
  char * resBin[world_size - 1];
  for (int i = 1; i < world_size; i++) {
    resBin[i - 1] = new char[resSize[i - 1]];
    MPI_Irecv(&resBin[i - 1][0], resSize[i - 1], MPI_BYTE, i, SENDAFFINCONT,
              MPI_COMM_WORLD, &requests[i - 1]);
  }
  MPI_Waitall(world_size - 1, requests, MPI_STATUS_IGNORE);
  info->infoMsg("Master got all binary data");
  // Deserialize binary messages
  for (int i = 1; i < world_size; i++) {
    std::vector<std::pair<std::string, float>> results;
    deserialize(results, resBin[i - 1], resSize[i - 1]);
    bucketResults.push_back(results);
    free(resBin[i - 1]);
  }
  info->infoMsg("Master deserialized everything");
  std::cout << "Master got following results:" << std::endl;
  for (int i = 1; i < world_size; i++) {
    std::cout << "From Worker #" << i << ":" << std::endl;
    for (auto j : bucketResults.at(i - 1)) {
      std::cout << j.first << ": " << j.second << std::endl;
    }
  }
  // Add the results to map and return FASTA sequences of added results
  std::vector<std::string> returnVal;
  for (int i = 1; i < world_size; i++) {
    for (unsigned int j = 0; j < bucketResults.at(i - 1).size(); j++) {
      std::string path = bucketResults.at(i - 1).at(j).first;
      size_t lastSlash = path.find_last_of("/");
      std::string prePath = path.substr(0, lastSlash);
      size_t secondToLastSlash = prePath.find_last_of("/");
      std::string fasta = prePath.substr(secondToLastSlash + 1,
                                         prePath.size() - secondToLastSlash);
      returnVal.push_back(fasta);
      float aff = bucketResults.at(i - 1).at(j).second;
      std::get<2>(internalMap[fasta]) = aff;
    }
  }
  return returnVal;
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
      throw PoolManagerException("Could not create directory for PDB file",
                                 FASTASEQ);
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
  return addElementsFromFiles(newFiles, world_size);
}

std::vector<std::string> PoolMGR::addElementsFromFASTAs(
                                    std::vector<std::string> &fastas,
                                    int world_size) {
  // Prepare files
  info->infoMsg("Fastas: " + std::to_string(fastas.size()));
  std::vector<std::string> newFiles;
  for (auto i : fastas) {
    if (internalMap.count(i) == 0) {
      internalMap[i] = std::make_tuple("", "", 10.0f, 0);
      genPDB(i);
      newFiles.push_back(std::get<1>(internalMap[i]));
    }
  }
  return addElementsFromFiles(newFiles, world_size);
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
