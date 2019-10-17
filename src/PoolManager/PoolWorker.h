/* Copyright 2019 iGEM Team Freiburg 2019
 *
 * Receives vector of FILES to perform Docking and MD on,
 * using OpenMP calculates affinities on available threads, returns
 * vector of <file, affinity> pairs
*/
#ifndef SRC_POOLMANAGER_POOLWORKER_H_
#define SRC_POOLMANAGER_POOLWORKER_H_
#include <omp.h>
#include <mpi.h>
#include <vector>
#include <string>
#include <utility>
#include "PoolManager.h"
#include "../Serialization/Serialization.h"
#include "../GMXInstance/GMXInstance.h"
#include "../Communication.h"
#include "../inih/INIReader.h"
#include "../VinaInstance/VinaInstance.h"
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

#endif  //  SRC_POOLMANAGER_POOLWORKER_H_
