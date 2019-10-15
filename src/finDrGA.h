/* Copyright 2019 Fabian Krause
 *
 * finDrGA Main class
 *
 * Tries to find a optimal ligand for given receptor(s) or multiple
 * conformations of receptor(s) using genetic algorithm strategies
*/
#ifndef SRC_FINDRGA_H_
#define SRC_FINDRGA_H_
#include <mpi.h>
#include <sys/stat.h>
#include <fstream>
#include "lib/GenAlgInst.h"
#include "finDrGAGenome.h"
#include "finDrGAFitnessFunc.h"
#include "PoolManager/PoolManager.h"
#include "inih/INIReader.h"
#include "cxxopts/cxxopts.hpp"
#endif  // SRC_FINDRGA_H_
