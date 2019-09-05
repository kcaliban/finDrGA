#include "PoolManager/PoolManager.h"
#include <sys/stat.h>
#include <limits.h>
#include <unistd.h>
#include <omp.h>

#include <stdlib.h>
#include <stdio.h>

std::string getexepath()
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  return std::string( result, (count > 0) ? count : 0 );
}

int main(int argc, char *argv[])
{
  const char * vinaPath = "/home/fk/autodock_vina_1_1_2_linux_x86/bin";
  const char * pythonShPath = "/home/fk/MGLTools-1.5.6/bin/pythonsh";
  const char * mgltoolstilitiesPath = "/home/fk/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24";
  const char * receptor = "/home/fk/Documents/iGEM/GenAlg/VinaGA/WorkingDir/AFAFA.pdb";
  const char * workDir = "/home/fk/Documents/iGEM/GenAlg/VinaGA/WorkingDir";

  PoolMGR poolMGR = PoolMGR(workDir, vinaPath, pythonShPath,
                            mgltoolstilitiesPath, receptor, 4, 5);

  std::vector<std::string> seqs = {"AFAFA", "FAFAF"};
  #pragma omp parallel
  #pragma omp for
  for (unsigned long i = 0; i < seqs.size(); i++) {
    // std::cout << "Adding sequence: " << seqs.at(i) << std::endl;
    poolMGR.addElement(seqs.at(i));
  }

  poolMGR.printSeqAff();
  return 0;
}
