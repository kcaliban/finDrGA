/* Copyright Fabian Krause 2019
 *
 * Receives vector of FILES to perform Docking and MD on,
 * using OpenMP calculates affinities on available processors, returns
 * vector of FILES, Affinity pairs
 *
 */
#include "PoolSlave.h"

std::string stripDir(std::string str) {
  size_t pos = str.find_last_of("/");
  return str.substr(0, pos);
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

int main (int argc, char **argv) {
  int world_size, world_rank;
  // Initialize the MPI environment
  MPI_Init(&argc , &argv);

  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Read in config
  while (42) {
    // Wait to receive a vector of files
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
    std::vector<std::pair<std::string, float>> results;
    #pragma omp parallel
    #pragma omp for
    for (unsigned int j = 0; j < FILES.size(); j++) {
    }


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
