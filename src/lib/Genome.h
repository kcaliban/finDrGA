/* Copyright 2019 iGEM Team Freiburg 2019
 *
 * Defines how individuals of GenoType are recombined and mutated in GA
*/
#ifndef SRC_LIB_GENOME_H_
#define SRC_LIB_GENOME_H_
#include <vector>
template <typename GenoType>
class Genome {
 public:
    // Crossover Function
    virtual GenoType crossOver(GenoType &, GenoType &, ...) = 0;
    // Mutation
    virtual GenoType mutate(GenoType &, ...) = 0;
};

#endif  // SRC_LIB_GENOME_H_
