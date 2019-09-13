/* Copyright 2019 Fabian Krause */
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
