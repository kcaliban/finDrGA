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
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Create topology using force field
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" pdb2gmx -f ");
  cmd.append(workDir);
  cmd.append("/clean.pdb");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/processed.gro");
  cmd.append(" -p ");
  cmd.append(workDir);
  cmd.append("/topol.top");
  cmd.append(" -i ");
  cmd.append(workDir);
  cmd.append("/posre.itp");
  cmd.append(" -water ");
  cmd.append(water);
  cmd.append(" -ff ");
  cmd.append(forcefield);
  cmd.append(" -ignh");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not generate topology for MD!\n"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Define the bounding box
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" editconf");
  cmd.append(" -f ");
  cmd.append(workDir);
  cmd.append("/processed.gro");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/newbox.gro");
  cmd.append(" -c ");
  cmd.append(" -d ");
  cmd.append(std::to_string(boxsize));
  cmd.append(" -bt ");
  cmd.append(bt);
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not define bounding box for MD!\n"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Solvate
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" solvate");
  cmd.append(" -cp ");
  cmd.append(workDir);
  cmd.append("/newbox.gro");
  cmd.append(" -cs ");
  cmd.append("spc216.gro");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/solv.gro");
  cmd.append(" -p ");
  cmd.append(workDir);
  cmd.append("/topol.top"); // Where will this file be?
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not solvate for MD!"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Add ions
  // Step one
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" grompp");
  cmd.append(" -f ");
  cmd.append(mdpPath);
  cmd.append("/ions.mdp");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/solv.gro");
  cmd.append(" -p ");
  cmd.append(workDir);
  cmd.append("/topol.top"); // Where will this file be?
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/ions.tpr");
  cmd.append(" -po ");
  cmd.append(workDir);
  cmd.append("/mdout.mdp");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not perform step one of ion adding"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Step two
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" genion");
  cmd.append(" -s ");
  cmd.append(workDir);
  cmd.append("/ions.tpr");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/solv_ions.gro");
  cmd.append(" -p ");
  cmd.append(workDir);
  cmd.append("/topol.top"); // Where will this file be?
  cmd.append(" -pname ");
  cmd.append("NA");
  cmd.append(" -nname ");
  cmd.append("CL");
  cmd.append(" -neutral ");
  cmd.append("<<eof\n13\neof"); // group SOL, might have to change to 16 depending
                                // on gromacs version
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not perform step two of ion adding"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Energy minimization
  // Prepare
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" grompp");
  cmd.append(" -f ");
  cmd.append(mdpPath);
  cmd.append("/minim.mdp");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/solv_ions.gro");
  cmd.append(" -p ");
  cmd.append(workDir);
  cmd.append("/topol.top");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/em.tpr");
  cmd.append(" -po ");
  cmd.append(workDir);
  cmd.append("/mdout.mdp");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not prepare energy minimzation"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Run MD for enery minimization
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" mdrun");
  cmd.append(" -gcom 2");
  cmd.append(" -s ");
  cmd.append(workDir);
  cmd.append("/em.tpr");
  cmd.append(" -deffnm em");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/em.gro");
  cmd.append(" -e ");
  cmd.append(workDir);
  cmd.append("/em.edr");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/em.trr");
  cmd.append(" -g ");
  cmd.append(workDir);
  cmd.append("/em.log");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not do energy minimzation"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Temperature Equilibrium
  // Preparation
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" grompp");
  cmd.append(" -f ");
  cmd.append(mdpPath);
  cmd.append("/nvt.mdp");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/em.gro");
  cmd.append(" -r ");
  cmd.append(workDir);
  cmd.append("/em.gro");
  cmd.append(" -p ");
  cmd.append(workDir);
  cmd.append("/topol.top");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/nvt.tpr");
  cmd.append(" -po ");
  cmd.append(workDir);
  cmd.append("/mdout.mdp");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not prepare establishing of equilibrium"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Run MD for equilibrium
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" mdrun");
  cmd.append(" -gcom 2");
  cmd.append(" -deffnm nvt");
  cmd.append(" -s ");
  cmd.append(workDir);
  cmd.append("/nvt.tpr");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/nvt.gro");
  cmd.append(" -e ");
  cmd.append(workDir);
  cmd.append("/nvt.edr");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/nvt.trr");
  cmd.append(" -cpo ");
  cmd.append(workDir);
  cmd.append("/nvt.cpt");
  cmd.append(" -g ");
  cmd.append(workDir);
  cmd.append("/nvt.log");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not establish equilibrim"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Pressure Equilibrium
  // Preparation
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" grompp");
  cmd.append(" -f ");
  cmd.append(mdpPath);
  cmd.append("/npt.mdp");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/nvt.gro");
  cmd.append(" -r ");
  cmd.append(workDir);
  cmd.append("/nvt.gro");
  cmd.append(" -t ");
  cmd.append(workDir);
  cmd.append("/nvt.cpt");
  cmd.append(" -p ");
  cmd.append(workDir);
  cmd.append("/topol.top");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/npt.tpr");
  cmd.append(" -po ");
  cmd.append(workDir);
  cmd.append("/mdout.mdp");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not prepare establishing of equilibrium"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Run MD for equilibrium
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" mdrun");
  cmd.append(" -gcom 2");
  cmd.append(" -deffnm npt");
  cmd.append(" -s ");
  cmd.append(workDir);
  cmd.append("/npt.tpr");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/npt.gro");
  cmd.append(" -e ");
  cmd.append(workDir);
  cmd.append("/npt.edr");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/npt.trr");
  cmd.append(" -g ");
  cmd.append(workDir);
  cmd.append("/npt.log");
  cmd.append(" -cpo ");
  cmd.append(workDir);
  cmd.append("/npt.cpt");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not establish equilibrim"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Final preparation
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" grompp");
  cmd.append(" -f ");
  cmd.append(mdpPath);
  cmd.append("/md.mdp");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/npt.gro");
  cmd.append(" -t ");
  cmd.append(workDir);
  cmd.append("/npt.cpt");
  cmd.append(" -p ");
  cmd.append(workDir);
  cmd.append("/topol.top");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/md_0_1.tpr");
  cmd.append(" -po ");
  cmd.append(workDir);
  cmd.append("/mdout.mdp");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Could not prepare MD tpr file"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
}

void GMXInstance::runMD() {
  /*
  // Run MD
  std::string cmd;
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" mdrun");
  cmd.append(" -deffnm md_0_1");
  cmd.append(" -s ");
  cmd.append(workDir);
  cmd.append("/md_0_1.tpr");
  cmd.append(" -c ");
  cmd.append(workDir);
  cmd.append("/md_0_1.gro");
  cmd.append(" -e ");
  cmd.append(workDir);
  cmd.append("/md_0_1.edr");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/md_0_1.trr");
  cmd.append(" -g ");
  cmd.append(workDir);
  cmd.append("/md_0_1.log");
  cmd.append(" -cpo ");
  cmd.append(workDir);
  cmd.append("/md_0_1.cpt");
  cmd.append(" -x ");
  cmd.append(workDir);
  cmd.append("/md_0_1.xtc");
  int success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Error running the MD!"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  */
  std::string cmd;
  int success;
  // Generate .pdb file
  // Step one
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" trjconv");
  cmd.append(" -s ");
  cmd.append(workDir);
  cmd.append("/md_0_1.tpr");
  cmd.append(" -f ");
  cmd.append(workDir);
  cmd.append("/md_0_1.xtc");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/md_0_1_noPBC.xtc");
  cmd.append(" -pbc mol ");
  cmd.append(" -center ");
  cmd.append("<<eof\n1\n0\neof");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Error trying to generate PDB from MD!"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Step two
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" trjconv");
  cmd.append(" -s ");
  cmd.append(workDir);
  cmd.append("/md_0_1.tpr");
  cmd.append(" -f ");
  cmd.append(workDir);
  cmd.append("/md_0_1_noPBC.xtc");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/MD.pdb");
  cmd.append(" <<eof\n1\neof");
  success = system(cmd.c_str());
  if (success == -1) {
    std::cout << "Error trying to generate PDB from MD!"
              << "Ligand: " << ligand
              << std::endl;
    exit(-1);
  }
  cmd.clear();
}

std::string GMXInstance::clusteredMD() {
  return "";
}
