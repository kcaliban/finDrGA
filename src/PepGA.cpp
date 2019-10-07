/* Copyright 2019 Fabian Krause */
#include "PepGA.h"

void check(const std::string p) {
  struct stat st;
  if (stat(p.c_str(), &st) != 0) {
    std::cout << "File or path \"" + p + "\" does not exist\n"
                   "Check your config.ini" << std::endl;
    exit(-1);
  }
}

void checkExecutable(const std::string e, const std::string progr) {
  std::string command;
  command.append(e);
  if (progr == "pymol") {
    command.append(" -kcq");
  } else if (progr == "vina") {
    command.append(" --help");
  } else if (progr == "pythonsh") {
    command.append(" -h");
  }
  command.append(" 2>/dev/null 1>&2");
  int success = system(command.c_str());
  if (success != 0) {
    std::cout << "Executable \"" + e + "\" does not exist or returns an error\n"
                   "Check your config.ini" << std::endl;
    exit(-1);
  }
}

std::string genToStr(std::vector<std::string> gen, PoolMGR * poolmgr) {
  std::string returnStr;
  returnStr.append("[");
  for (auto g : gen) {
    returnStr.append(g);
    returnStr.append(": ");
    returnStr.append(std::to_string(poolmgr->getAffinity(g)));
    returnStr.append(", ");
  }
  returnStr.pop_back(); returnStr.pop_back();
  returnStr.append("]");

  return returnStr;
}

// Get PDBs from initialpop folder
std::vector<std::string> getInitialPop(std::string dir) {
  std::string command;
  command.append("ls ");
  command.append(dir);
  command.append(" | cat");
  std::string output;
  FILE * lsOutputStream = popen(command.c_str(), "r");
  char buf[1024];
  while (fgets(buf, 1024, lsOutputStream)) {
    output += buf;
  }
  pclose(lsOutputStream);

  std::vector<std::string> result;
  std::string line;
  std::stringstream outputStream(output);
  while (std::getline(outputStream, line, '\n')) {
    result.push_back(line);
  }
  return result;
}

// Get a random sample of PDB filenames from specified folder
std::vector<std::string> getRandomSample(std::string dir, int amount) {
  std::string command;
  command.append("ls ");
  command.append(dir);
  command.append(" | shuf -n ");
  command.append(std::to_string(amount));
  std::string output;
  FILE * lsOutputStream = popen(command.c_str(), "r");
  char buf[1024];
  while (fgets(buf, 1024, lsOutputStream)) {
    output += buf;
  }
  pclose(lsOutputStream);
  // Create and return vector
  std::vector<std::string> result;
  std::string line;
  std::stringstream outputStream(output);
  while (std::getline(outputStream, line, '\n')) {
    result.push_back(line);
  }
  return result;
}

// Get a single PDB
std::string getRandomPDB(std::string dir) {
  return getRandomSample(dir, 1).at(0);
}

// Get receptor filenames
std::vector<std::string> getReceptorsM(std::string dir, bool prep = false) {
  // Read all pdb files in a directory
  std::string command;
  command.append("ls ");
  command.append(dir);
  command.append(" | cat ");
  if (!prep) {
    command.append("| grep -v .pdbqt | grep -v conf");
  } else {
    command.append("| grep -v conf");
  }
  std::string output;
  FILE * lsOutputStream = popen(command.c_str(), "r");
  char buf[1024];
  while (fgets(buf, 1024, lsOutputStream)) {
    output += buf;
  }
  pclose(lsOutputStream);
  // Return vector of every filename
  std::vector<std::string> result;
  std::string line;
  std::stringstream outputStream(output);
  while (std::getline(outputStream, line, '\n')) {
    result.push_back(line);
  }
  return result;
}

void prepareConfig(std::string receptor) {
  // Generate conf file for docking

  // Read receptor file
  std::ifstream t(receptor);
  t.seekg(0, std::ios::end);
  size_t size = t.tellg();
  std::string buffer(size, ' ');
  t.seekg(0);
  t.read(&buffer[0], size);
  t.close();

  // std::cout << "Read receptor file!" << std::endl;

  // For max and min no sorting is required
  float xmin = std::numeric_limits<float>::infinity();
  float ymin = std::numeric_limits<float>::infinity();
  float zmin = std::numeric_limits<float>::infinity();
  float xmax = - std::numeric_limits<float>::infinity();
  float ymax = - std::numeric_limits<float>::infinity();
  float zmax = - std::numeric_limits<float>::infinity();

  // Read line per line, set min and max accordingly
  std::string line;
  std::stringstream receptorStream(buffer);
  while (std::getline(receptorStream, line, '\n')) {
    if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
      float x = stof(line.substr(31, 8));
      xmin = (x < xmin) ? x : xmin;
      xmax = (x > xmax) ? x : xmax;
      float y = stof(line.substr(39, 8));
      ymin = (y < ymin) ? y : ymin;
      ymax = (y > ymax) ? y : ymax;
      float z = stof(line.substr(47, 8));
      zmin = (z < zmin) ? z : zmin;
      zmax = (z > zmax) ? z : zmax;
    }
  }

  float sizex = xmax - xmin;
  float sizey = ymax - ymin;
  float sizez = zmax - zmin;

  // float offsetx = xmin + sizex / 2;
  float offsetx = xmin + sizex / 2.0;
  float offsety = ymin + sizey / 2.0;
  float offsetz = zmin + sizez / 2.0;

  std::string outfile = receptor + "_conf";

  std::ofstream confFile;
  confFile.open(outfile.c_str(), std::ios::trunc);
  if (!confFile) {
    std::cout << "Could not open config file!" << std::endl;
    exit(-1);
  }
  confFile << "center_x = " << offsetx << std::endl;
  confFile << "center_y = " << offsety << std::endl;
  confFile << "center_z = " << offsetz << std::endl;
  confFile << std::endl;
  confFile << "size_x = " << sizex + 30 << std::endl;
  confFile << "size_y = " << sizey + 30 << std::endl;
  confFile << "size_z = " << sizez + 30 << std::endl;
  confFile.close();
}

void preparePDBQT(std::string receptor,
                   std::string pythonShPath,
                   std::string mgltoolstilitiesPath) {
  // Generate a PDBQT
  std::string command;
  command.append(pythonShPath);
  command.append(" ");
  command.append(mgltoolstilitiesPath);
  command.append("/prepare_receptor4.py -r ");
  command.append(receptor);
  command.append(" -A bonds_hydrogens -U nphs -o ");
  command.append(receptor);
  command.append("qt");
  int success = system(command.c_str());
  if (success != 0) {
    throw VinaException("Could not generate pdbqt file for receptor", receptor);
  }
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cout << "Wrong number of arguments!" << std::endl;
    std::cout << "Usage: PepGA [number of generations]"
                 " [prob. of random point mutation for each individual]"
                 " [percentage of top individuals to copy each gen.]"
                 " [size of initial pop., taken from initialpdbs"
                 " and randompdbs if initialpdbs does not suffice]"
      << std::endl;
    return 1;
  }
  unsigned int noPop = atoi(argv[1]);
  float mutateProb = atof(argv[2]);
  float genCpy = atof(argv[3]);
  unsigned int gen = (argc == 5) ? atoi(argv[4]) : 10;  // def. size 10
  /* Initialize OpenMPI */
  int world_size, world_rank;
  // Initialize the MPI environment
  MPI_Init(&argc , &argv);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  // Print out information about main process
  /**************/
  /* Prepare random engine */
  std::random_device rd;
  std::mt19937 mt(rd());
  /**************/
  /* Read config */
  INIReader reader("config.ini");
  if (reader.ParseError() != 0) {
        std::cout << "Can't load 'config.ini'\n"
                     "Check if it exists in the same dir as PepGA";
        return 1;
  }
  // Executables
  std::string vinaPath = reader.Get("paths", "vina", "vina");
  checkExecutable(vinaPath, "vina");
  std::string pythonShPath = reader.Get("paths", "pythonsh", "pythonsh");
  checkExecutable(pythonShPath, "pythonsh");
  std::string pymolPath = reader.Get("paths", "pymol", "pymol");
  checkExecutable(pymolPath, "pymol");
  std::string gromacsPath = reader.Get("paths", "gromacs", "gmx");
  checkExecutable(gromacsPath, "");
  // Required by Vina/preparation for Vina
  std::string mgltoolstilitiesPath = reader.Get("paths", "MGLToolsUtilities",
                                                "");
  check(mgltoolstilitiesPath);
  std::string workDir = reader.Get("paths", "workingDir", "");
  check(workDir);
  std::string receptorsPath = reader.Get("paths", "receptors", "");
  check(receptorsPath);
  bool receptorsPrep = reader.GetBoolean("paths", "receptorsprep", false);
  int exhaustiveness = reader.GetInteger("VINA", "exhaustiveness", 1);
  int energy_range = reader.GetInteger("VINA", "energy_range", 5);
  // Required by gromacs
  std::string settings = reader.Get("paths", "settings", "");
  check(settings);
  std::string forcefieldPath = reader.Get("paths", "forcefieldpath", "");
  check(forcefieldPath);
  std::string forcefield = reader.Get("GROMACS", "forcefield", "");
  std::string water = reader.Get("GROMACS", "water", "");
  std::string boundingboxtype = reader.Get("GROMACS", "bt", "");
  float boxsize = reader.GetReal("GROMACS", "boxsize", 1.0);
  float clustercutoff = reader.GetReal("GROMACS", "clustercutoff", 0.12);
  // Path to PDBs for first generation
  std::string initialpdbs = reader.Get("paths", "initialpdbs", "");
  bool initialpdbsMD = reader.GetBoolean("paths", "initialpdbsisMD", false);
  // Path to PDB files to take random sample from
  std::string randompdbs = reader.Get("paths", "randompdbs", "");
  if (!initialpdbs.empty()) {check(initialpdbs);}
  if (!randompdbs.empty()) {check(randompdbs);}
  /**************/
  Info info(true, true, workDir + "/" + "PepLOG");
  /* Print info about master node */
  info.infoMsg("Master has rank " + std::to_string(world_rank) + "(should be 0)");
  /* Get receptors */
  std::vector<std::string> receptors;
  std::vector<std::string> receptorfiles = getReceptorsM(receptorsPath,
                                                        receptorsPrep);
  for (auto s : receptorfiles) {
    receptors.push_back(receptorsPath + "/" + s);
  }
  if (!receptorsPrep) {
    for (auto s : receptors) {
      prepareConfig(s);
      try {
        preparePDBQT(s, pythonShPath, mgltoolstilitiesPath);
      } catch (std::exception& e) {
        info.errorMsg(e.what(), true);
      }
    }
  }
  /**************/
  /* Generate ligands */
  // Initialization of key objects required
  GenAlgInst<std::string, PepGenome, PepFitnessFunc> inst(&mt);
  PoolMGR poolmgr(workDir.c_str(), vinaPath.c_str(), pythonShPath.c_str(),
                  mgltoolstilitiesPath.c_str(), pymolPath.c_str(),
                  receptors,
                  exhaustiveness, energy_range, gromacsPath.c_str(),
                  settings.c_str(), forcefield.c_str(), forcefieldPath.c_str(),
                  water.c_str(), boundingboxtype.c_str(), boxsize,
                  clustercutoff, &info);
  PepFitnessFunc fitnessFunc(&poolmgr);
  PepGenome vinaGenome(&mt);
  // Initial pdbs
  std::vector<std::string> startingSequences;
  info.infoMsg("Gathering the initial population...");
  std::vector<std::string> initPopulation;
  if (initialpdbs != "") {
    info.infoMsg("Adding peptides from initialpdbs & "
                 "randompdbs to the gene pool...");
    std::vector<std::string> initPDB = getInitialPop(initialpdbs);
    for (auto i : initPDB) {
      initPopulation.push_back(initialpdbs + "/" + i);
    }
  }
  if (randompdbs != "" && gen - initPopulation.size() > 0) {
    std::vector<std::string> randomSample = getRandomSample(randompdbs,
                                                 gen - initPopulation.size());
    for (auto i : randomSample) {
      initPopulation.push_back(randompdbs + "/" + i);
    }
  }
  info.errorMsg("Population size in PoolMGR: " + std::to_string(initPopulation.size()), false);
  startingSequences = poolmgr.addElementsFromPDBs(initPopulation, world_size);
  /**************/
  /* GA */
  std::vector<std::string> curGen = startingSequences;
  info.infoMsg("POPULATION SIZE: " + std::to_string(curGen.size()));
  for (unsigned int i = 0; i < noPop; i++) {
    // Output to log file
    std::string output = "Generation: ";
    output.append(std::to_string(i));
    output.append("\nIndividuals:\n");
    output.append(genToStr(curGen, &poolmgr));
    output.append("\nItems in Pool Manager:\n");
    output.append(poolmgr.toStr());
    info.infoMsg(output);
    // Get new generation
    curGen = inst.nextGen(vinaGenome, fitnessFunc, curGen, mutateProb, genCpy);

    // Remove duplicates for adding
    std::vector<std::string> curGenDistinct;
    curGenDistinct = curGen;
    // Sort and then remove consecutive duplicates
    std::sort(curGenDistinct.begin(), curGenDistinct.end());
    curGenDistinct.erase(std::unique(curGenDistinct.begin(),
                                     curGenDistinct.end()),
                          curGenDistinct.end());
    // Add the new elements
    try {
      poolmgr.addElementsFromFASTAs(curGenDistinct, world_size);
    } catch (std::exception& e) {
      info.errorMsg(e.what(), true);
    }
  }
  /**************/
  MPI_Finalize();
  return 0;
}
