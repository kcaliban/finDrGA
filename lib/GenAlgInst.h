#ifndef GENMA
#define GENMA
#include "Genome.h"
#include "FitnessFunction.h"
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <fstream>
template <typename GenoType, typename Genome, typename FitnessFunction>
class GenAlgInst
{
  public:
    void simulate(Genome genome, FitnessFunction fitnessfunc,
                                     std::vector<GenoType> genotype, int n,
                                     float mutateProb, bool debug=true,
                                     bool entropy=true,
                                     char * entropyFile="entropy") {
      std::vector<GenoType> curGen = genotype;
      for (int i = 0; i < n; i++) {
        if (debug) {
          std::cout << "Generation: " << (i + 1) << std::endl;
        }
        // Will the old gen be deleted automatically?
        curGen = nextGen(genome, fitnessfunc, curGen, mutateProb,
                          debug, entropy, entropyFile);
      }
    };
  private:
    int calculateEntropy(std::vector<GenoType> genotypes) {
      unsigned int entropy = 0;
      // Naive O(n^2) Worst Case implementation
      for (unsigned int i = 0; i < genotypes.size(); i++) {
        bool unique = true;
        for (unsigned int j = i + 1; j < genotypes.size(); j++) {
          if (genotypes[i] == genotypes[j]) {
            unique = false;
            break;
          }
        }
        if (unique) {
          entropy++;
        }
      }
      return entropy;
    }

    std::vector<GenoType> nextGen(Genome genome, FitnessFunction fitnessfunc,
                                    std::vector<GenoType> genotypes,
                                    float mutateProb, bool debug=true,
                                    bool entropy=true, char * entropyFile="entropy") {
      std::vector<GenoType> newGen;
      // Required for selection and recombination
      std::vector<float> fitnesses;
      if (debug) {
        std::cout << "\tCalculating fitnesses..." << std::endl;
      }
      for (GenoType genotype : genotypes) {
        float fitness = fitnessfunc.calculateFitness(genotype);
        fitnesses.push_back(fitness);
      }
      // SELECTION (this should become "modular" aswell)
      //           (the way selection and recombination is done should be
      //            defined elsewhere)
      // Copy the top 20%
      // Requires them to be sorted; we need to keep the original order
      // of elements to associate them with genotypes vector (could also
      // use pairs)
      if (debug) {
        std::cout << "\tSelection..." << std::endl;
      }
      std::vector<size_t> sortedindices(fitnesses.size());
      std::iota(sortedindices.begin(), sortedindices.end(), 0);
      sort(sortedindices.begin(), sortedindices.end(),
          [fitnesses](size_t i1, size_t i2) {
            return fitnesses[i1] > fitnesses[i2];});
      int amount = (int) (0.2 * genotypes.size());
      for (int i = 0; i < amount; i++) {
        newGen.push_back(genotypes[sortedindices[i]]);
      }
      if (debug) {
        std::cout << "\tBest individual: " << genotypes[sortedindices[0]]
                  << std::endl;
      }
      // RECOMBINATION
      // Generate a discrete random distribution for fitnesses, i.e.
      // each element has prob. of its weight (fitness) divided by the sum
      // of all weights
      if (debug) {
        std::cout << "\tRecombination..." << std::endl;
      }
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator(seed);
      std::discrete_distribution<int> fitnessdistribution(fitnesses.begin(), fitnesses.end());
      // Pick two genotypes randomly, until we have a population as big as the initial
      while (amount < genotypes.size()) {
        if (debug) {
          std::cout << "\t\tPopulation size: " << amount << std::endl;
        }
        GenoType inda = genotypes[fitnessdistribution(generator)];
        GenoType indb = genotypes[fitnessdistribution(generator)];
        newGen.push_back(genome.crossOver(inda, indb));
        amount++;
      }
      // MUTATION
      // With probability mutateProb mutate every individual
      // Use a uniform distribution between 0 & 1, if random generated number
      // is <= mutateProb do a mutation
      if (debug) {
        std::cout << "\tMutation..." << std::endl;
      }
      std::uniform_real_distribution<float> uniformdistribution(0.0, 1.0);
      for (int i = 0; i < newGen.size(); i++) {
        if (uniformdistribution(generator) <= mutateProb) {
          GenoType mutatedGen = genome.mutate(newGen[i]);
          /*
          // Save old genotype for deletion
          GenoType old = newGen.at(i);
          */
          // Insert the new one
          newGen.at(i) = mutatedGen;
          // delete old;
        }
      }
      if (entropy) {
          std::ofstream outfile;

          outfile.open(entropyFile, std::ios::out | std::ios::app);
          outfile.close();
          outfile.open(entropyFile, std::ios_base::app);
          outfile << calculateEntropy(newGen) << std::endl;
      }

      // std::move?
      return newGen;
    }
};

#endif
