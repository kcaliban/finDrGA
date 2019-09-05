#include "VinaGA.h"

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout << "Wrong number of arguments!" << std::endl;
    std::cout << "Usage: VinaGA [Number of populations]"
      << std::endl;
    exit(-1);
  }
  const char * vinaPath = "/home/fk/autodock_vina_1_1_2_linux_x86/bin";
  const char * pythonShPath = "/home/fk/MGLTools-1.5.6/bin/pythonsh";
  const char * mgltoolstilitiesPath = "/home/fk/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24";
  const char * receptor = "/home/fk/Documents/iGEM/GenAlg/VinaGA/WorkingDir/AFAFA.pdb";
  const char * workDir = "/home/fk/Documents/iGEM/GenAlg/VinaGA/WorkingDir";

  GenAlgInst<std::string, VinaGenome, VinaFitnessFunc> inst;
  std::vector<std::string> startingSequences = {"AFG", "YFA", "AYE", "HGA",
      "BAY", "AAF", "AAG", "HHA", "HGG", "HAG", "AGF", "BYA", "YYA",
      "AYA", "AAY"};
  PoolMGR poolmgr(workDir, vinaPath, pythonShPath, mgltoolstilitiesPath,
                   receptor, 4, 5);

  // Initialize Pool manager for starting sequences
  #pragma omp parallel
  #pragma omp for
  for (unsigned long i = 0; i < startingSequences.size(); i++) {
    // std::cout << "Adding sequence: " << seqs.at(i) << std::endl;
    poolmgr.addElement(startingSequences.at(i));
  }

  VinaFitnessFunc fitnessFunc(&poolmgr);
  VinaGenome vinaGenome;

  inst.simulate(vinaGenome, fitnessFunc, startingSequences, 2000, 0.05);

  std::cout << "lol" << std::endl;
  poolmgr.printSeqAff();
  return 0;
}
