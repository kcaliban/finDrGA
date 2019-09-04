#include "GMXInstance.h"
#include <sys/stat.h>
#include <limits.h>
#include <unistd.h>
#include <omp.h>

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  const char * gromacsPath = "/usr/local/gromacs/bin";
  const char * ligand = "/home/fk/Documents/iGEM/GenAlg/VinaGA/GMXInstance/work/AFAFA.pdb";
  const char * workDir = "/home/fk/Documents/iGEM/GenAlg/VinaGA/GMXInstance/work";
  const char * mdpPath = "/home/fk/Documents/iGEM/GenAlg/VinaGA/GMXInstance/templ";
  const char * forcefield = "amber99sb-ildn-fme"; // how to do with correct directory?
  const char * water = "spce";
  const char * bt = "dodecahedron";
  float boxsize = 1.0;

  GMXInstance gmxInstance(ligand, gromacsPath, workDir, forcefield,
                          water, bt, boxsize, mdpPath);

  // gmxInstance.preparePDB();
  gmxInstance.runMD();
  return 0;
}
