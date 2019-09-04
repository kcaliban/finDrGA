#include "GMXInstance.h"

void GMXInstance::preparePDB() {
  // Clean file from crystal water
  std::string cmd;
  cmd.append("grep -v HOH ");
  cmd.append(ligand);
  cmd.append(" > ");
  cmd.append(workDir);
  cmd.append("/");
  cmd.append("clean.pdb");
  int success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not clean PDB file for MD!\n"
              << "Ligand name: " << ligand
              << std::endl;
    exit(-1);
  }
  // Create topology using force field
  
}
