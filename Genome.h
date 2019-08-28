#ifndef GENOM
#define GENOM
#include <vector>
// Singleton? Then no coevolution possible.
template <typename GenoType>
class Genome {
  public:
    // Crossover Function
    virtual GenoType * crossOver(GenoType, GenoType, ...) = 0;
    // Mutation
    virtual GenoType * mutate(GenoType, ...) = 0;
};

#endif
