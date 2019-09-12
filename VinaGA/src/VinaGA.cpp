#include "VinaGA.h"

std::string genToStr(std::vector<std::string> gen, PoolMGR poolmgr) {
  std::string returnStr;
  returnStr.append("[");
  for (auto g : gen) {
    returnStr.append(g);
    returnStr.append(": ");
    returnStr.append(std::to_string(poolmgr.getAffinity(g)));
    returnStr.append(", ");
  }
  returnStr.pop_back(); returnStr.pop_back();
  returnStr.append("]");

  return returnStr;
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

// Get receptor filenames
std::vector<std::string> getReceptors(std::string dir) {
  // Read all pdb files in a directory
  std::string command;
  command.append("ls ");
  command.append(dir);
  command.append(" | cat | grep -v .pdbqt");
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

int main(int argc, char *argv[])
{
  if (argc != 3 && argc != 4) {
    std::cout << "Wrong number of arguments!" << std::endl;
    std::cout << "Usage: VinaGA [Number of populations]"
                 " [prob. of random poinmutation]"
                 " [optional: size of initial gen (if initialpdbs specified)]"
      << std::endl;
    exit(-1);
  }

  INIReader reader("config.ini");
  if (reader.ParseError() != 0) {
        std::cout << "Can't load 'config.ini'\n";
        return 1;
  }

  // AutoDock VINA
  std::string vinaPath = reader.Get("paths", "vina", "vina");
  std::string pythonShPath = reader.Get("paths", "pythonsh", "pythonsh");
  std::string mgltoolstilitiesPath = reader.Get("paths", "MGLToolsUtilities", "");
  std::string pymolPath = reader.Get("paths", "pymol", "pymol");
  std::string workDir = reader.Get("paths", "workingDir", "");
  std::string receptorsPath = reader.Get("paths", "receptors", "");
  int exhaustiveness = reader.GetInteger("VINA", "exhaustiveness", 1);
  int energy_range = reader.GetInteger("VINA", "energy_range", 5);

  // GROMACS
  std::string gromacsPath = reader.Get("paths", "gmx", "gmx");
  std::string mdpPath = reader.Get("paths", "mdp", "");
  std::string forcefieldPath = reader.Get("paths", "forcefieldpath", "");
  std::string forcefield = reader.Get("GROMACS", "forcefield", "");
  std::string water = reader.Get("GROMACS", "water", "");
  std::string boundingboxtype = reader.Get("GROMACS", "bt", "");
  float boxsize = reader.GetReal("GROMACS", "bt", 1.0);
  float clustercutoff = reader.GetReal("GROMACS", "clustercutoff", 0.12);

  // Path to PDB files to take random sample from
  std::string initialpdbs = reader.Get("paths", "initialpdbs", "");
  // Standard number of randomly picked individuals: 10
  int gen = 10;
  if (argc != 4) {
    // Size of initial gen not specified
    initialpdbs = "";
  } else {
    gen = atoi(argv[3]);
  }

  // Get receptors
  std::vector<std::string> receptors;
  std::vector<std::string> receptorfiles = getReceptors(receptorsPath);
  for (auto s : receptorfiles) {
    receptors.push_back(receptorsPath + "/" + s);
  }

  GenAlgInst<std::string, VinaGenome, VinaFitnessFunc> inst;
  PoolMGR poolmgr(workDir.c_str(), vinaPath.c_str(), pythonShPath.c_str(),
                  mgltoolstilitiesPath.c_str(), pymolPath.c_str(),
                  receptors,
                  exhaustiveness, energy_range, gromacsPath.c_str(),
                  mdpPath.c_str(), forcefield.c_str(), forcefieldPath.c_str(),
                  water.c_str(), boundingboxtype.c_str(), boxsize,
                  clustercutoff);
  VinaFitnessFunc fitnessFunc(&poolmgr);
  VinaGenome vinaGenome;

  // Gather elements
  std::vector<std::string> startingSequences;
  if (initialpdbs == "") {
    // Just some random sequences for testing
    startingSequences = {
         "HLYE", "LAFY", "IAGY", "YHVL", "AHGG", "KPAG", "HAGF", "BYAH", "CYLA"};
  } else {
    // Gather random sample of files
    std::vector<std::string> files = getRandomSample(initialpdbs, gen);
    // Add each file to pool whilst adding their fasta sequences to initial population
    // for later MD and Docking
    for (auto filename : files) {
      std::string FASTA = poolmgr.addElementPDB(initialpdbs + "/" + filename);
      if (!FASTA.empty()) {
        startingSequences.push_back(FASTA);
      }
    }
  }

  std::vector<std::string> curGen = startingSequences;
  for (int i = 0; i < atoi(argv[1]); i++) {
    // Remove duplicates for adding
    std::vector<std::string> curGenDistinct;
    curGenDistinct = curGen;
    // Sort and then remove consecutive duplicates
    std::sort(curGenDistinct.begin(), curGenDistinct.end());
    curGenDistinct.erase(std::unique(curGenDistinct.begin(), curGenDistinct.end()),
                          curGenDistinct.end());
    // Add the new elements (PoolMGR only adds them if they don't exist already, does a new MD else)
    #pragma omp parallel
    #pragma omp for
    for (unsigned int i = 0; i < curGenDistinct.size(); i++) {
      poolmgr.addElement(curGenDistinct.at(i));
    }
    // CleanUp of not used strings
    // poolmgr.update(curGen);
    // poolmgr.cleanUp(3);
    // Output to log file
    std::string output = "Generation: ";
    output.append(std::to_string(i));
    output.append("\n");
    output.append(genToStr(curGen, poolmgr));
    output.append("\n");
    std::string outputFileDir = workDir;
    outputFileDir.append("/VINAGALOG");
    std::ofstream outputFileStream(outputFileDir,  std::ios::out | std::ios::app);
    outputFileStream << output;
    outputFileStream.close();
    // Output for best/entropy graph
    // Get new generation
    curGen = inst.nextGen(vinaGenome, fitnessFunc, curGen, atof(argv[2]));
    // Temporary debug output
    poolmgr.printSeqAff();
  }

  return 0;
}
