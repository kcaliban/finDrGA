/* Copyright Fabian Krause 2019
 *
 * Receives vector of FILES to perform Docking and MD on,
 * using OpenMP calculates affinities on available processors, returns
 * vector of <file, affinity> pairs
 *
 */
#include "PoolWorker.h"

std::string stripDir(std::string str) {
  size_t pos = str.find_last_of("/");
  return str.substr(0, pos);
}

void preparePDBQT(std::string ligand) {
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

float genDock(std::string file) {
  float affinity = 10;
  std::string fileCluster = stripDir(file) + "/topcluster.pdb";
  // Prepare ligand
  try {
    preparePDBQT(fileCluster);
  } catch (...) {
    throw;
  }
  // Do a docking for each receptor
  for (unsigned int i = 0; i < receptors.size(); i++) {
    VinaInstance vinaInstance(vinaPath.c_str(), receptors.at(i).c_str(),
                              fileCluster.c_str(),
                              info);
    float recaffinity = vinaInstance.calculateBindingAffinity(exhaustiveness,
                                                           energy_range);
    if (recaffinity < affinity) { affinity = recaffinity; }
  }

  return affinity;
}

void genMD(std::string file) {
  GMXInstance gmxInstance(file.c_str(),
                          gromacsPath.c_str(),
                          pymolPath.c_str(),
                          stripDir(file).c_str(),
                          forcefield.c_str(),
                          forcefieldPath.c_str(),
                          water.c_str(),
                          boundingboxtype.c_str(),
                          clustercutoff,
                          boxsize,
                          mdpPath.c_str(),
                          info);
  try {
    gmxInstance.preparePDB();
    gmxInstance.runMD();
    gmxInstance.clusterMD();
    gmxInstance.extractTopCluster();
  } catch (...) {
    throw;
  }
}

// Get receptor filenames
std::vector<std::string> getReceptors(std::string dir, bool prep = false) {
  // Read all pdb files in a directory
  std::string command;
  command.append("ls ");
  command.append(dir);
  command.append(" | cat ");
  if (!prep) {
    command.append("| grep -v .pdbqt | grep -v conf");
  } else {
    command.append("| grep -v conf");
  }
  std::string output;
  FILE * lsOutputStream = popen(command.c_str(), "r");
  char buf[1024];
  while (fgets(buf, 1024, lsOutputStream)) {
    output += buf;
  }
  pclose(lsOutputStream);
  // Return vector of every filename
  std::vector<std::string> result;
  std::string line;
  std::stringstream outputStream(output);
  while (std::getline(outputStream, line, '\n')) {
    result.push_back(line);
  }
  return result;
}

int main (int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(&argc , &argv);

  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Read in relevant parameters
  INIReader reader("config.ini");
  if (reader.ParseError() != 0) {
        std::cout << "Can't load 'config.ini'\n"
                     "Check if it exists in the same dir as PepGA";
        return 1;
  }
  pymolPath = reader.Get("paths", "pymol", "pymol");
  gromacsPath = reader.Get("paths", "gromacs", "gmx");
  forcefield = reader.Get("GROMACS", "forcefield", "");
  forcefieldPath = reader.Get("paths", "forcefieldpath", "");
  water = reader.Get("GROMACS", "water", "");
  boundingboxtype = reader.Get("GROMACS", "bt", "");
  boxsize = reader.GetReal("GROMACS", "boxsize", 1.0);
  clustercutoff = reader.GetReal("GROMACS", "clustercutoff", 0.12);
  mdpPath = reader.Get("paths", "settings", "");
  exhaustiveness = reader.GetInteger("VINA", "exhaustiveness", 1);
  energy_range = reader.GetInteger("VINA", "energy_range", 5);
  vinaPath = reader.Get("paths", "vina", "vina");
  pythonShPath = reader.Get("paths", "pythonsh", "pythonsh");
  mgltoolstilitiesPath = reader.Get("paths", "MGLToolsUtilities",
                                                "");
  std::string receptorsPath = reader.Get("paths", "receptors", "");
  info = new Info(false, true, ""); // Console output

  std::vector<std::string> oldRec = getReceptors(receptorsPath, false);
  for (auto i : oldRec) {
    receptors.push_back(receptorsPath + "/" + i);
  }
  while (42) {
    // Wait to receive a vector of files
    info->infoMsg("Worker #" + std::to_string(world_rank)
                             + " waiting for a job...");
    std::vector<std::string> FILES;
    unsigned int FILESSize;

    MPI_Recv(&FILESSize, 1, MPI_INT, 0, SENDFILESSIZE,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    char * tmp = new char[FILESSize];
    MPI_Recv(&tmp[0], FILESSize, MPI_BYTE, 0, SENDFILESCONT,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    deserialize(FILES, tmp, FILESSize);
    free(tmp);

    // Do the right thing
    info->infoMsg("Worker #" + std::to_string(world_rank) + " got a job!");
    std::vector<std::pair<std::string, float>> results;
    #pragma omp parallel
    #pragma omp for
    for (unsigned int j = 0; j < FILES.size(); j++) {
      // Do MD
      try {
        genMD(FILES.at(j));
      } catch (...) {
        info->errorMsg("MD for " + FILES.at(j) +
                       " failed, skipping...", false);
      }
      // Do Docking
      float aff = genDock(FILES.at(j));
      // Add result
      try {
        results.push_back(std::make_pair(FILES.at(j), aff));
      } catch (...) {
        info->errorMsg("Docking for " + FILES.at(j) +
                       " failed, skipping...", false);
      }
    }
    info->infoMsg("Worker #" + std::to_string(world_rank)
                             + " is sending back the results now!");

    // Send back the results
    unsigned int resultsSize;
    tmp = serialize(results, &resultsSize);
    MPI_Send(&resultsSize, 1, MPI_INT, 0, SENDAFFINSIZE, MPI_COMM_WORLD);
    MPI_Send(&tmp[0], resultsSize, MPI_BYTE, 0, SENDAFFINCONT, MPI_COMM_WORLD);
    free(tmp);
  }

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
