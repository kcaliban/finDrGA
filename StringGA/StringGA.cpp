#include "StringGA.h"
int main(int argc, char *argv[])
{
  if (argc < 5) {
    std::cout << "Wrong number of arguments!" << std::endl;
    std::cout << "Usage: StringGA [Goal string] [Initial population size]"
      " [Number of generations] [Prob. of mutation] {optional: Name of entropy file}" << std::endl;
    exit(-1);
  }
  bool entropy = false;
  char * entropyFile = {};
  if (argc == 6) {
    entropy = true;
    entropyFile = argv[5];
  }
  GenAlgInst<std::string,
             StringGenome,
             StringFitnessFunc>
    inst;
  int popSize = atoi(argv[2]);
  std::string goalString = argv[1];
  // Generate random strings
  std::string alphabet =
                        " ,"
                        "0123456789"
                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                        "abcdefghijklmnopqrstuvwxyz";
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_int_distribution<int> distribution (0, alphabet.size() - 1);
  std::vector<std::string> initGen;
  for (int i = 0; i < popSize; i++) {
    std::string newString;
    for (size_t j = 0; j < goalString.size(); j++) {
      newString.push_back(alphabet[distribution(generator)]);
    }
    initGen.push_back(newString);
  }
  std::cout << "Initial generation: [";
  for (auto gen : initGen) {
    std::cout << gen << ",";
  }
  std::cout << "]" << std::endl;
  // Create our fitness function
  StringFitnessFunc fitnessFunc = StringFitnessFunc(goalString);
  // "Create" a Genome (this is ugly)
  StringGenome genome = StringGenome(alphabet);
  std::cout << strtod(argv[3], nullptr) << std::endl;
  // Run the simulation
  inst.simulate(genome, fitnessFunc, initGen, atoi(argv[3]), strtod(argv[4], nullptr),
      true, entropy, entropyFile);
}
