#include "VinaGA.h"

int main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cout << "Wrong number of arguments!" << std::endl;
    std::cout << "Usage: VinaGA [Number of populations]"
                 " [prob. of random poinmutation]"
      << std::endl;
    exit(-1);
  }

  INIReader reader("config.ini");
  if (reader.ParseError() != 0) {
        std::cout << "Can't load 'config.ini'\n";
        return 1;
  }

  std::string vinaPath = reader.Get("paths", "vina", "");
  std::string pythonShPath = reader.Get("paths", "pythonsh", "");
  std::string mgltoolstilitiesPath = reader.Get("paths", "MGLToolsUtilities", "");
  std::string workDir = reader.Get("paths", "workingDir", "");
  std::string receptor = reader.Get("paths", "receptor", "");
  int exhaustiveness = reader.GetInteger("options", "exhaustiveness", 1);
  int energy_range = reader.GetInteger("options", "energy_range", 5);

  GenAlgInst<std::string, VinaGenome, VinaFitnessFunc> inst;
  std::vector<std::string> startingSequences = {"AFG", "YFA", "AYE", "HGA"};
     // "BAY", "AAF", "AAG", "HHA", "HGG", "HAG", "AGF", "BYA", "YYA",
     // "AYA", "AAY"};

  PoolMGR poolmgr(workDir.c_str(), vinaPath.c_str(), pythonShPath.c_str(),
                  mgltoolstilitiesPath.c_str(), receptor.c_str(),
                  exhaustiveness, energy_range);

  VinaFitnessFunc fitnessFunc(&poolmgr);
  VinaGenome vinaGenome;

  std::vector<std::string> curGen = startingSequences;
  for (int i = 0; i < atoi(argv[1]); i++) {
    // Add the new elements (PoolMGR only adds them if they don't exist already)
    #pragma omp parallel
    #pragma omp for
    for (unsigned int i = 0; i < curGen.size(); i++) {
      poolmgr.addElement(curGen.at(i));
    }
    // CleanUp of not used strings
    poolmgr.update(curGen);
    poolmgr.cleanUp(3);
    // Get new generation
    curGen = inst.nextGen(vinaGenome, fitnessFunc, curGen, atof(argv[2]));
  }

  poolmgr.printSeqAff();
  return 0;
}
