/* Copyright 2019 Fabian Krause
 *
 * PepGA Main class
 *
 * Tries to find a optimal ligand for given receptor(s) or multiple
 * conformations of receptor(s) using genetic algorithm strategies
*/
#ifndef SRC_PEPGA_H_
#define SRC_PEPGA_H_
#include <mpi.h>
#include <fstream>
#include <sys/stat.h>
#include "lib/GenAlgInst.h"
#include "PepGenome.h"
#include "PepFitnessFunc.h"
#include "PoolManager/PoolManager.h"
#include "inih/INIReader.h"
#endif  // SRC_PEPGA_H_
