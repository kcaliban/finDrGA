/* Copyright 2019 Fabian Krause
 *
 * Dvelopr Main class
 *
 * Tries to find a optimal ligand for given receptor(s) or multiple
 * conformations of receptor(s) using genetic algorithm strategies
*/
#ifndef SRC_DVELOPR_H_
#define SRC_DVELOPR_H_
#include <mpi.h>
#include <sys/stat.h>
#include <fstream>
#include "lib/GenAlgInst.h"
#include "DveloprGenome.h"
#include "DveloprFitnessFunc.h"
#include "PoolManager/PoolManager.h"
#include "inih/INIReader.h"
#include "cxxopts/cxxopts.hpp"
#endif  // SRC_DVELOPR_H_
