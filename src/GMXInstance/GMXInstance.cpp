#include "GMXInstance.h"

std::string GMXInstance::logStr() {
  if (log) {
    return " >> " + workDir + "/GMXINSTLOG" + " 2>&1";
  }
  return " > /dev/null 2>&1";
}

void GMXInstance::debugPrint(const char * str) {
  if (debug) {
    std::cout << "\033[1;32mINFO (MD, " << ligand << "): " << str << "\033[0m" << std::endl;
  }
}

void GMXInstance::preparePDB() {
  // Export path to forcefield
  debugPrint("Setting env forcefield value...");
  std::string command;
  int success = setenv("GMXLIB", forcefieldPath.c_str(), 1);
  if (success != 0) {
    throw (GMXException("Could not set GMXLib Path", ""));
  }
  command.clear();
  // Clean file from crystal water
  debugPrint("Cleaning ligand from crystal water...");
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
  // Create topology using force field
  debugPrint("Creating topology...");
  command.append(gromacsPath);
  command.append(" pdb2gmx -f ");
  command.append(workDir);
  command.append("/clean.pdb");
  command.append(" -o ");
  command.append(workDir);
  command.append("/processed.gro");
  command.append(" -p ");
  command.append(workDir);
  command.append("/topol.top");
  command.append(" -i ");
  command.append(workDir);
  command.append("/posre.itp");
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
  debugPrint("Defining the bounding box...");
  command.append(gromacsPath);
  command.append(" editconf");
  command.append(" -f ");
  command.append(workDir);
  command.append("/processed.gro");
  command.append(" -o ");
  command.append(workDir);
  command.append("/newbox.gro");
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
  debugPrint("Solvating...");
  command.append(gromacsPath);
  command.append(" solvate");
  command.append(" -cp ");
  command.append(workDir);
  command.append("/newbox.gro");
  command.append(" -cs ");
  command.append("spc216.gro");
  command.append(" -o ");
  command.append(workDir);
  command.append("/solv.gro");
  command.append(" -p ");
  command.append(workDir);
  command.append("/topol.top"); // Where will this file be?
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not solvate for MD", ligand);
  }
  command.clear();
  // Add ions
  debugPrint("Adding ions...");
  // Step one
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/ions.mdp");
  command.append(" -c ");
  command.append(workDir);
  command.append("/solv.gro");
  command.append(" -p ");
  command.append(workDir);
  command.append("/topol.top");
  command.append(" -o ");
  command.append(workDir);
  command.append("/ions.tpr");
  command.append(" -po ");
  command.append(workDir);
  command.append("/mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not ionize for MD (1)", ligand);
  }
  command.clear();
  // Step two
  command.append(gromacsPath);
  command.append(" genion");
  command.append(" -s ");
  command.append(workDir);
  command.append("/ions.tpr");
  command.append(" -o ");
  command.append(workDir);
  command.append("/solv_ions.gro");
  command.append(" -p ");
  command.append(workDir);
  command.append("/topol.top"); // Where will this file be?
  command.append(" -pname ");
  command.append("NA");
  command.append(" -nname ");
  command.append("CL");
  command.append(" -neutral ");
  command.append(logStr());
  command.append(" ");
  command.append("<<eof\n13\neof"); // group SOL, might have to change to 16 depending
                                // on gromacs version
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not ionize for MD (2)", ligand);
  }
  command.clear();
  // Energy minimization
  debugPrint("Minimzing energy...");
  // Prepare
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/minim.mdp");
  command.append(" -c ");
  command.append(workDir);
  command.append("/solv_ions.gro");
  command.append(" -p ");
  command.append(workDir);
  command.append("/topol.top");
  command.append(" -o ");
  command.append(workDir);
  command.append("/em.tpr");
  command.append(" -po ");
  command.append(workDir);
  command.append("/mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare energy minimzation", ligand);
  }
  command.clear();
  // Run MD for enery minimization
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -gcom 2");
  command.append(" -s ");
  command.append(workDir);
  command.append("/em.tpr");
  command.append(" -deffnm em");
  command.append(" -c ");
  command.append(workDir);
  command.append("/em.gro");
  command.append(" -e ");
  command.append(workDir);
  command.append("/em.edr");
  command.append(" -o ");
  command.append(workDir);
  command.append("/em.trr");
  command.append(" -g ");
  command.append(workDir);
  command.append("/em.log");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not do energy minimzation", ligand);
  }
  command.clear();
  // Temperature Equilibrium
  debugPrint("Equilibriating temperature...");
  // Preparation
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/nvt.mdp");
  command.append(" -c ");
  command.append(workDir);
  command.append("/em.gro");
  command.append(" -r ");
  command.append(workDir);
  command.append("/em.gro");
  command.append(" -p ");
  command.append(workDir);
  command.append("/topol.top");
  command.append(" -o ");
  command.append(workDir);
  command.append("/nvt.tpr");
  command.append(" -po ");
  command.append(workDir);
  command.append("/mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare establishing of equilibrium", ligand);
  }
  command.clear();
  // Run MD for equilibrium
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -gcom 2");
  command.append(" -deffnm nvt");
  command.append(" -s ");
  command.append(workDir);
  command.append("/nvt.tpr");
  command.append(" -c ");
  command.append(workDir);
  command.append("/nvt.gro");
  command.append(" -e ");
  command.append(workDir);
  command.append("/nvt.edr");
  command.append(" -o ");
  command.append(workDir);
  command.append("/nvt.trr");
  command.append(" -cpo ");
  command.append(workDir);
  command.append("/nvt.cpt");
  command.append(" -g ");
  command.append(workDir);
  command.append("/nvt.log");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not establish equilibrim", ligand);
  }
  command.clear();
  // Pressure Equilibrium
  debugPrint("Equilibriating pressure...");
  // Preparation
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/npt.mdp");
  command.append(" -c ");
  command.append(workDir);
  command.append("/nvt.gro");
  command.append(" -r ");
  command.append(workDir);
  command.append("/nvt.gro");
  command.append(" -t ");
  command.append(workDir);
  command.append("/nvt.cpt");
  command.append(" -p ");
  command.append(workDir);
  command.append("/topol.top");
  command.append(" -o ");
  command.append(workDir);
  command.append("/npt.tpr");
  command.append(" -po ");
  command.append(workDir);
  command.append("/mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare establishing of equilibrium", ligand);
  }
  command.clear();
  // Run MD for equilibrium
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -gcom 2");
  command.append(" -deffnm npt");
  command.append(" -s ");
  command.append(workDir);
  command.append("/npt.tpr");
  command.append(" -c ");
  command.append(workDir);
  command.append("/npt.gro");
  command.append(" -e ");
  command.append(workDir);
  command.append("/npt.edr");
  command.append(" -o ");
  command.append(workDir);
  command.append("/npt.trr");
  command.append(" -g ");
  command.append(workDir);
  command.append("/npt.log");
  command.append(" -cpo ");
  command.append(workDir);
  command.append("/npt.cpt");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not establish equilibrim", ligand);
  }
  command.clear();
  // Final preparation
  debugPrint("Final preparation for MD...");
  command.append(gromacsPath);
  command.append(" grompp");
  command.append(" -f ");
  command.append(mdpPath);
  command.append("/md.mdp");
  command.append(" -c ");
  command.append(workDir);
  command.append("/npt.gro");
  command.append(" -t ");
  command.append(workDir);
  command.append("/npt.cpt");
  command.append(" -p ");
  command.append(workDir);
  command.append("/topol.top");
  command.append(" -o ");
  command.append(workDir);
  command.append("/md_0_1.tpr");
  command.append(" -po ");
  command.append(workDir);
  command.append("/mdout.mdp");
  command.append(logStr());
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not prepare MD tpr file", ligand);
  }
  command.clear();
}

void GMXInstance::runMD() {
  // Run MD
  debugPrint("Running the MD...");
  std::string command;
  command.append(gromacsPath);
  command.append(" mdrun");
  command.append(" -deffnm md_0_1");
  command.append(" -s ");
  command.append(workDir);
  command.append("/md_0_1.tpr");
  command.append(" -c ");
  command.append(workDir);
  command.append("/md_0_1.gro");
  command.append(" -e ");
  command.append(workDir);
  command.append("/md_0_1.edr");
  command.append(" -o ");
  command.append(workDir);
  command.append("/md_0_1.trr");
  command.append(" -g ");
  command.append(workDir);
  command.append("/md_0_1.log");
  command.append(" -cpo ");
  command.append(workDir);
  command.append("/md_0_1.cpt");
  command.append(" -x ");
  command.append(workDir);
  command.append("/md_0_1.xtc");
  command.append(logStr());
  int success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not run the MD", ligand);
  }
  command.clear();
  debugPrint("MD successful!");
  // Generate .pdb file
  // Step one
  debugPrint("Generating PDB file...");
  command.append(gromacsPath);
  command.append(" trjconv");
  command.append(" -s ");
  command.append(workDir);
  command.append("/md_0_1.tpr");
  command.append(" -f ");
  command.append(workDir);
  command.append("/md_0_1.xtc");
  command.append(" -o ");
  command.append(workDir);
  command.append("/md_0_1_noPBC.xtc");
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
  command.append(gromacsPath);
  command.append(" trjconv");
  command.append(" -s ");
  command.append(workDir);
  command.append("/md_0_1.tpr");
  command.append(" -f ");
  command.append(workDir);
  command.append("/md_0_1_noPBC.xtc");
  command.append(" -o ");
  command.append(workDir);
  command.append("/MD.pdb");
  command.append(logStr());
  command.append(" ");
  command.append(" <<eof\n1\neof");
  success = system(command.c_str());
  if (success != 0) {
    throw GMXException("Could not generate PDB from MD", ligand);
  }
  command.clear();
}

void GMXInstance::clusteredMD() {
  std::string command;
  debugPrint("Clustering PDB...");
  command.append(gromacsPath);
  command.append(" cluster");
  command.append(" -f ");
  command.append(workDir);
  command.append("/MD.pdb");
  command.append(" -o ");
  command.append(workDir);
  command.append("/rmsd-clust.xpm");
  command.append(" -g ");
  command.append(workDir);
  command.append("/clust.log");
  command.append(" -s ");
  command.append(workDir);
  command.append("/md_0_1.tpr");
  command.append(" -dist ");
  command.append(workDir);
  command.append("/rmsd-dist.xvg");
  command.append(" -cl ");
  command.append(workDir);
  command.append("/clusters.pdb");
  command.append(" -sz ");
  command.append(workDir);
  command.append("/clust-size.xvg");
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
  debugPrint("Extracting the top cluster...");
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
  std::string line;
  std::stringstream clustSizeFileStream(buffer);
  std::regex clusterSizeRegEx("^\\s*([0-9]+)\\s*([0-9]+)\\s*$");// Regex for number
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
    throw GMXException("Could not extract top cluster from clustered MD", ligand);
  }
}
