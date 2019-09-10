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

void GMXInstance::errorPrint(const char * str) {
  std::cout << "\033[1;31mERROR (MD): " << str << std::endl;
  std::cout << "Ligand: " << ligand << "\033[0m" << std::endl;

  exit(-1);
}

void GMXInstance::preparePDB() {
  // Export path to forcefield
  debugPrint("Setting env forcefield value...");
  std::string cmd;
  int success = setenv("GMXLIB", forcefieldPath.c_str(), 1);
  if (success != 0) {
    std::cout << "Error trying to set GMXLib Path!"
              << std::endl;
    exit(-1);
  }
  cmd.clear();
  // Clean file from crystal water
  debugPrint("Cleaning ligand from crystal water...");
  cmd.append("grep -v HOH ");
  cmd.append(ligand);
  cmd.append(" > ");
  cmd.append(workDir);
  cmd.append("/");
  cmd.append("clean.pdb");
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not clean PDB file for MD!");
  }
  cmd.clear();
  // Create topology using force field
  debugPrint("Creating topology...");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not generate topology for MD!");
  }
  cmd.clear();
  // Define the bounding box
  debugPrint("Defining the bounding box...");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not define bounding box for MD!");
  }
  cmd.clear();
  // Solvate
  debugPrint("Solvating...");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not solvate for MD!");
  }
  cmd.clear();
  // Add ions
  debugPrint("Adding ions...");
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
  cmd.append("/topol.top");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/ions.tpr");
  cmd.append(" -po ");
  cmd.append(workDir);
  cmd.append("/mdout.mdp");
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not perform step one of ion adding");
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
  cmd.append(logStr());
  cmd.append(" ");
  cmd.append("<<eof\n13\neof"); // group SOL, might have to change to 16 depending
                                // on gromacs version
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not perform step two of ion adding");
  }
  cmd.clear();
  // Energy minimization
  debugPrint("Minimzing energy...");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not prepare energy minimzation");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not do energy minimzation");
  }
  cmd.clear();
  // Temperature Equilibrium
  debugPrint("Equilibriating temperature...");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not prepare establishing of equilibrium");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not establish equilibrim");
  }
  cmd.clear();
  // Pressure Equilibrium
  debugPrint("Equilibriating pressure...");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not prepare establishing of equilibrium");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not establish equilibrim");
  }
  cmd.clear();
  // Final preparation
  debugPrint("Final preparation for MD...");
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
  cmd.append(logStr());
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not prepare MD tpr file");
  }
  cmd.clear();
}

void GMXInstance::runMD() {
  // Run MD
  debugPrint("Running the MD...");
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
  cmd.append(logStr());
  int success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Error running the MD!");
  }
  cmd.clear();
  debugPrint("MD successful!");
  // Generate .pdb file
  // Step one
  debugPrint("Generating PDB file...");
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
  cmd.append(logStr());
  cmd.append(" ");
  cmd.append("<<eof\n1\n0\neof");
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Error trying to generate PDB from MD!");
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
  cmd.append(logStr());
  cmd.append(" ");
  cmd.append(" <<eof\n1\neof");
  success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Error trying to generate PDB from MD!");
  }
  cmd.clear();
}

void GMXInstance::clusteredMD() {
  std::string cmd;
  debugPrint("Clustering PDB...");
  cmd.append(gromacsPath);
  cmd.append("/gmx");
  cmd.append(" cluster");
  cmd.append(" -f ");
  cmd.append(workDir);
  cmd.append("/MD.pdb");
  cmd.append(" -o ");
  cmd.append(workDir);
  cmd.append("/rmsd-clust.xpm");
  cmd.append(" -g ");
  cmd.append(workDir);
  cmd.append("/clust.log");
  cmd.append(" -s ");
  cmd.append(workDir);
  cmd.append("/md_0_1.tpr");
  cmd.append(" -dist ");
  cmd.append(workDir);
  cmd.append("/rmsd-dist.xvg");
  cmd.append(" -cl ");
  cmd.append(workDir);
  cmd.append("/clusters.pdb");
  cmd.append(" -sz ");
  cmd.append(workDir);
  cmd.append("/clust-size.xvg");
  cmd.append(" -cutoff ");
  cmd.append(std::to_string(clustercutoff));
  cmd.append(" ");
  cmd.append(logStr());
  cmd.append(" <<eof\n1\n1\neof");
  int success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Error clustering the MD!");
  }
  cmd.clear();
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
  std::string cmd;
  cmd.append(pymolPath);
  cmd.append(" -kcQ -d \"load ");
  cmd.append(workDir);
  cmd.append("/clusters.pdb;");
  cmd.append(" save ");
  cmd.append(workDir);
  cmd.append("/topcluster.pdb, state=");
  cmd.append(std::to_string(maxNo));
  cmd.append("\"");
  cmd.append(logStr());
  int success = system(cmd.c_str());
  if (success != 0) {
    errorPrint("Could not extract top cluster from clustered MD!");
  }
}
