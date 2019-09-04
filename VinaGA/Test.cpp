#include "PoolManager.h"
#include <sys/stat.h>
#include <limits.h>
#include <unistd.h>

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
  // const char * receptor = "AFAFA_D.pdb";
  const char * receptor = "/home/fk/Documents/iGEM/GenAlg/VinaGA/AFAFA.pdb";

  const char * workDir = "/home/fk/Documents/iGEM/GenAlg/VinaGA/WorkingDir";

  PoolMGR poolMGR = PoolMGR(workDir, vinaPath, pythonShPath,
                            mgltoolstilitiesPath, receptor, 1, 5);
  poolMGR.addElement("AFAFA");
  poolMGR.addElement("FAFAF");
  poolMGR.printSeqAff();
  return 0;
}
