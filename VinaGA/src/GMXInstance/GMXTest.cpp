#include "GMXInstance.h"
#include <sys/stat.h>
#include <limits.h>
#include <unistd.h>

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  const char * gromacsPath = "/usr/local/gromacs/bin";
  const char * ligand = "/home/fk/Documents/iGEM/GenAlg/VinaGA/src/GMXInstance/work/AFAFA.pdb";
  const char * workDir = "/home/fk/Documents/iGEM/GenAlg/VinaGA/src/GMXInstance/work";
  const char * mdpPath = "/home/fk/Documents/iGEM/GenAlg/VinaGA/src/GMXInstance/templ";
  const char * forcefield = "amber99sb-ildn-fme"; // how to do with correct directory?
  const char * forcefieldPath = "/home/fk/Documents/iGEM/GenAlg/VinaGA/src/GMXInstance/forcefields";
  const char * water = "spce";
  const char * bt = "dodecahedron";
  float boxsize = 1.0;

  GMXInstance gmxInstance(ligand, gromacsPath, workDir, forcefield, forcefieldPath,
                          water, bt, boxsize, mdpPath, true, true);

  gmxInstance.preparePDB();
  gmxInstance.runMD();
  return 0;
}
