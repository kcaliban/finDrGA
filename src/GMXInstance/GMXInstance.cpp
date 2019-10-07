/* Copyright 2019 Fabian Krause */
#include "GMXInstance.h"

std::string GMXInstance::logStr() {
  return " >> " + workDir + "/GMXINSTLOG" + " 2>&1";
}

void GMXInstance::energyMinim() {
  // Export path to forcefield
  info->infoMsg("(GMX, " + ligand + ") Setting env forcefield value...");
  std::string command;
  int success = setenv("GMXLIB", forcefieldPath.c_str(), 1);
  if (success != 0) {
    throw(GMXException("Could not set GMXLib Path", ""));
  }
  command.clear();
  // Prepare for GROMACS
  info->infoMsg("(GMX, " + ligand + ") Preparing cleansed PDB for GROMACS...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" pdb2gmx -f ");
  command.append(ligand);
  command.append(" -o ");
  command.append("processed.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -i ");
  command.append("posre.itp");
  command.append(" -water ");
  command.append(water);
  command.append(" -ff ");
  command.append(forcefield);
  command.append(" -ignh");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not generate topology for MD", ligand, "TOP");
  }
  command.clear();
  // Define the bounding box
  info->infoMsg("(GMX, " + ligand + ") Defining the bounding box...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" editconf");
  command.append(" -f ");
  command.append("processed.gro");
  command.append(" -o ");
  command.append("newbox.gro");
  command.append(" -c ");
  command.append(" -d ");
  command.append(std::to_string(boxsize));
  command.append(" -bt ");
  command.append(bt);
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not define bounding box for MD", ligand);
  }
  command.clear();
  // Solvate
  info->infoMsg("(GMX, " + ligand + ") Solvating...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" solvate");
  command.append(" -cp ");
  command.append("newbox.gro");
  command.append(" -cs ");
  command.append("spc216.gro");
  command.append(" -o ");
  command.append("solv.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not solvate for MD", ligand);
  }
  command.clear();
  // Add ions
  info->infoMsg("(GMX, " + ligand + ") Adding ions...");
  // Step one
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/ions.mdp");
  command.append(" -c ");
  command.append("solv.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -o ");
  command.append("ions.tpr");
  command.append(" -po ");
  command.append("mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not ionize for MD (1)", ligand);
  }
  command.clear();
  // Step two
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" genion");
  command.append(" -s ");
  command.append("ions.tpr");
  command.append(" -o ");
  command.append("solv_ions.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -pname ");
  command.append("NA");
  command.append(" -nname ");
  command.append("CL");
  command.append(" -neutral ");
  command.append(logStr());
  command.append(" ");
  command.append("<<eof\n13\neof");  // group SOL, might have to change
                                     // to 16 depending on gromacs version
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not ionize for MD (2)", ligand);
  }
  command.clear();
  // Energy minimization
  info->infoMsg("(GMX, " + ligand + ") Minimzing energy...");
  // Prepare
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/minim.mdp");
  command.append(" -c ");
  command.append("solv_ions.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -o ");
  command.append("em.tpr");
  command.append(" -po ");
  command.append("mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare energy minimzation", ligand);
  }
  command.clear();
  // Run MD for energy minimization
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -nt 1");
  command.append(" -s ");
  command.append("em.tpr");
  command.append(" -deffnm em");
  command.append(" -c ");
  command.append("em.gro");
  command.append(" -e ");
  command.append("em.edr");
  command.append(" -o ");
  command.append("em.trr");
  command.append(" -g ");
  command.append("em.log");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not do energy minimzation", ligand);
  }
  command.clear();
  // Centering
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" trjconv");
  command.append(" -s ");
  command.append("em.tpr");
  command.append(" -f ");
  command.append("em.trr");
  command.append(" -o ");
  command.append("noPBC.xtc");
  command.append(" -pbc ");
  command.append("mol");
  command.append(" -center ");
  command.append(logStr());
  command.append(" <<eof\n1\n0\neof");
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not do energy minimzation", ligand);
  }
  command.clear();
  // gmx2pdb
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" trjconv");
  command.append(" -s ");
  command.append("em.tpr");
  command.append(" -f ");
  command.append("noPBC.xtc");
  command.append(" -o ");
  command.append("em.pdb");
  command.append(logStr());
  command.append(" <<eof\n1\neof");
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not do energy minimzation", ligand);
  }
  command.clear();
}

void GMXInstance::preparePDB() {
  // Export path to forcefield
  info->infoMsg("(GMX, " + ligand + ") Setting env forcefield value...");
  std::string command;
  int success = setenv("GMXLIB", forcefieldPath.c_str(), 1);
  if (success != 0) {
    throw(GMXException("Could not set GMXLib Path", ""));
  }
  command.clear();
  // Clean file from crystal water
  info->infoMsg("(GMX, " + ligand + ") Cleaning ligand from crystal water...");
  command.append("grep -v HOH ");
  command.append(ligand);
  command.append(" > ");
  command.append(workDir);
  command.append("/");
  command.append("clean.pdb");
  success = system(command.c_str());
  if (success != 0) {
    throw(GMXException("Could not clean PDB file for MD", ligand));
  }
  command.clear();
  // Prepare for GROMACS
  info->infoMsg("(GMX, " + ligand + ") Preparing cleansed PDB for GROMACS...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" pdb2gmx -f ");
  command.append("clean.pdb");
  command.append(" -o ");
  command.append("processed.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -i ");
  command.append("posre.itp");
  command.append(" -water ");
  command.append(water);
  command.append(" -ff ");
  command.append(forcefield);
  command.append(" -ignh");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not generate topology for MD", ligand, "TOP");
  }
  command.clear();
  // Define the bounding box
  info->infoMsg("(GMX, " + ligand + ") Defining the bounding box...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" editconf");
  command.append(" -f ");
  command.append("processed.gro");
  command.append(" -o ");
  command.append("newbox.gro");
  command.append(" -c ");
  command.append(" -d ");
  command.append(std::to_string(boxsize));
  command.append(" -bt ");
  command.append(bt);
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not define bounding box for MD", ligand);
  }
  command.clear();
  // Solvate
  info->infoMsg("(GMX, " + ligand + ") Solvating...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" solvate");
  command.append(" -cp ");
  command.append("newbox.gro");
  command.append(" -cs ");
  command.append("spc216.gro");
  command.append(" -o ");
  command.append("solv.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not solvate for MD", ligand);
  }
  command.clear();
  // Add ions
  info->infoMsg("(GMX, " + ligand + ") Adding ions...");
  // Step one
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/ions.mdp");
  command.append(" -c ");
  command.append("solv.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -o ");
  command.append("ions.tpr");
  command.append(" -po ");
  command.append("mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not ionize for MD (1)", ligand);
  }
  command.clear();
  // Step two
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" genion");
  command.append(" -s ");
  command.append("ions.tpr");
  command.append(" -o ");
  command.append("solv_ions.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -pname ");
  command.append("NA");
  command.append(" -nname ");
  command.append("CL");
  command.append(" -neutral ");
  command.append(logStr());
  command.append(" ");
  command.append("<<eof\n13\neof");  // group SOL, might have to change
                                     // to 16 depending on gromacs version
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not ionize for MD (2)", ligand);
  }
  command.clear();
  // Energy minimization
  info->infoMsg("(GMX, " + ligand + ") Minimzing energy...");
  // Prepare
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/minim.mdp");
  command.append(" -c ");
  command.append("solv_ions.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -o ");
  command.append("em.tpr");
  command.append(" -po ");
  command.append("mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare energy minimzation", ligand);
  }
  command.clear();
  // Run MD for enery minimization
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -nt 1");
  command.append(" -s ");
  command.append("em.tpr");
  command.append(" -deffnm em");
  command.append(" -c ");
  command.append("em.gro");
  command.append(" -e ");
  command.append("em.edr");
  command.append(" -o ");
  command.append("em.trr");
  command.append(" -g ");
  command.append("em.log");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not do energy minimzation", ligand);
  }
  command.clear();
  // Temperature Equilibrium
  info->infoMsg("(GMX, " + ligand + ") Equilibriating temperature...");
  // Preparation
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/nvt.mdp");
  command.append(" -c ");
  command.append("em.gro");
  command.append(" -r ");
  command.append("em.gro");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -o ");
  command.append("nvt.tpr");
  command.append(" -po ");
  command.append("mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare establishing of equilibrium", ligand);
  }
  command.clear();
  // Run MD for equilibrium
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -nt 1");
  command.append(" -deffnm nvt");
  command.append(" -s ");
  command.append("nvt.tpr");
  command.append(" -c ");
  command.append("nvt.gro");
  command.append(" -e ");
  command.append("nvt.edr");
  command.append(" -o ");
  command.append("nvt.trr");
  command.append(" -cpo ");
  command.append("nvt.cpt");
  command.append(" -g ");
  command.append("nvt.log");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not establish equilibrim", ligand);
  }
  command.clear();
  // Pressure Equilibrium
  info->infoMsg("(GMX, " + ligand + ") Equilibriating pressure...");
  // Preparation
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/npt.mdp");
  command.append(" -c ");
  command.append("nvt.gro");
  command.append(" -r ");
  command.append("nvt.gro");
  command.append(" -t ");
  command.append("nvt.cpt");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -o ");
  command.append("npt.tpr");
  command.append(" -po ");
  command.append("mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare establishing of equilibrium", ligand);
  }
  command.clear();
  // Run MD for equilibrium
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -nt 1");
  command.append(" -deffnm npt");
  command.append(" -s ");
  command.append("npt.tpr");
  command.append(" -c ");
  command.append("npt.gro");
  command.append(" -e ");
  command.append("npt.edr");
  command.append(" -o ");
  command.append("npt.trr");
  command.append(" -g ");
  command.append("npt.log");
  command.append(" -cpo ");
  command.append("npt.cpt");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not establish equilibrim", ligand);
  }
  command.clear();
  // Final preparation
  info->infoMsg("(GMX, " + ligand + ") Final preparation for MD...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/md.mdp");
  command.append(" -c ");
  command.append("npt.gro");
  command.append(" -t ");
  command.append("npt.cpt");
  command.append(" -p ");
  command.append("topol.top");
  command.append(" -o ");
  command.append("md_0_1.tpr");
  command.append(" -po ");
  command.append("mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare MD tpr file", ligand);
  }
  command.clear();
}

void GMXInstance::runMD() {
  // Run MD
  info->infoMsg("(GMX, " + ligand + ") Running the MD...");
  std::string command;
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -nt 1");
  command.append(" -deffnm md_0_1");
  command.append(" -s ");
  command.append("md_0_1.tpr");
  command.append(" -c ");
  command.append("md_0_1.gro");
  command.append(" -e ");
  command.append("md_0_1.edr");
  command.append(" -o ");
  command.append("md_0_1.trr");
  command.append(" -g ");
  command.append("md_0_1.log");
  command.append(" -cpo ");
  command.append("md_0_1.cpt");
  command.append(" -x ");
  command.append("md_0_1.xtc");
  command.append(logStr());
  int success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not run the MD", ligand);
  }
  command.clear();
  info->infoMsg("(GMX, " + ligand + ") MD successful!");
  // Generate .pdb file
  // Step one
  info->infoMsg("(GMX, " + ligand + ") Generating PDB file...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" trjconv");
  command.append(" -s ");
  command.append("md_0_1.tpr");
  command.append(" -f ");
  command.append("md_0_1.xtc");
  command.append(" -o ");
  command.append("md_0_1_noPBC.xtc");
  command.append(" -pbc mol ");
  command.append(" -center ");
  command.append(logStr());
  command.append(" ");
  command.append("<<eof\n1\n0\neof");
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not generate PDB from MD", ligand);
  }
  command.clear();
  // Step two
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" trjconv");
  command.append(" -s ");
  command.append("md_0_1.tpr");
  command.append(" -f ");
  command.append("md_0_1_noPBC.xtc");
  command.append(" -o ");
  command.append("MD.pdb");
  command.append(logStr());
  command.append(" ");
  command.append(" <<eof\n1\neof");
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not generate PDB from MD", ligand);
  }
  command.clear();
}

void GMXInstance::clusterMD() {
  std::string command;
  info->infoMsg("(GMX, " + ligand + ") Clustering PDB...");
  command.append("cd ");
  command.append(workDir);
  command.append("; ");
  command.append(gromacsPath);
  command.append(" cluster");
  command.append(" -f ");
  command.append("MD.pdb");
  command.append(" -o ");
  command.append("rmsd-clust.xpm");
  command.append(" -g ");
  command.append("clust.log");
  command.append(" -s ");
  command.append("md_0_1.tpr");
  command.append(" -dist ");
  command.append("rmsd-dist.xvg");
  command.append(" -cl ");
  command.append("clusters.pdb");
  command.append(" -sz ");
  command.append("clust-size.xvg");
  command.append(" -cutoff ");
  command.append(std::to_string(clustercutoff));
  command.append(" ");
  command.append(logStr());
  command.append(" <<eof\n1\n1\neof");
  int success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not cluster the MD", ligand);
  }
  command.clear();
}

void GMXInstance::extractTopCluster() {
  info->infoMsg("(GMX, " + ligand + ") Extracting the top cluster...");
  std::string clustSizeFile;
  clustSizeFile.append(workDir);
  clustSizeFile.append("/clust-size.xvg");

  // Read clust-size.xvg
  std::ifstream t(clustSizeFile);
  t.seekg(0, std::ios::end);
  size_t size = t.tellg();
  std::string buffer(size, ' ');
  t.seekg(0);
  t.read(&buffer[0], size);
  t.close();

  // Get the biggest cluster
  int maxNo;
  int max = - std::numeric_limits<int>::infinity();
  try {
    std::string line;
    std::stringstream clustSizeFileStream(buffer);
    std::regex clusterSizeRegEx("^\\s*([0-9]+)\\s*([0-9]+)\\s*$");
    std::smatch lineClusterMatch;
    while (std::getline(clustSizeFileStream, line, '\n')) {
      std::regex_search(line, lineClusterMatch, clusterSizeRegEx);
      if (lineClusterMatch.empty()) {continue;}
      int cur = stoi(lineClusterMatch.str(2));
      int curNo = stoi(lineClusterMatch.str(1));
      if (cur > max) {
        max = cur;
        maxNo = curNo;
      }
    }
  } catch (std::regex_error &e) {
    info->errorMsg(e.what(), false);
    info->errorMsg("Error in extracting the top cluster."
                   "Check clust-size.xvg and readjust your settings.", false);
    info->errorMsg("Contents of clust-size.xvg:", false);
    info->errorMsg(buffer, true);
  }


  // Extract biggest cluster
  std::string command;
  command.append(pymolPath);
  command.append(" -kcQ -d \"load ");
  command.append(workDir);
  command.append("/clusters.pdb;");
  command.append(" save ");
  command.append(workDir);
  command.append("/topcluster.pdb, state=");
  command.append(std::to_string(maxNo));
  command.append("\"");
  command.append(logStr());
  int success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not extract top cluster from clustered MD",
                       ligand);
  }
}
