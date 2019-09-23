#ifndef SRC_POOLMANAGER_POOLSLV_H_
#define SRC_POOLMANAGER_POOLSLV_H_
#include <omp.h>
#include <mpi.h>
#include <string>
#include <utility>
#include "PoolManager.h"
#include "../Serialization/Serialization.h"
#include "../GMXInstance/GMXInstance.h"
#include "../Communication.h"
#include "../inih/INIReader.h"
std::string pymolPath;
std::string gromacsPath;
std::string forcefield;
std::string forcefieldPath;
std::string water;
std::string boundingboxtype;
std::string mdpPath;
std::string vinaPath;
std::string pythonShPath;
std::string mgltoolstilitiesPath;
std::vector<std::string> receptors;

float boxsize;
float clustercutoff;
int exhaustiveness;
int energy_range;

Info * info;
int world_size, world_rank;

#endif  /* ifndef SRC_POOLMANAGER_POOLSLV_H_ */
