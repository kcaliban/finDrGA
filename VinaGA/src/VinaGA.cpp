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
  command.append(" | cat | grep -v .pdbqt | grep -v conf");
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
  // Size x,y,z are simply the max minus the min

  /*
  std::cout << "Trying to read receptor file: "
            << receptor << std::endl;
            */
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
    if (line.substr(0, 4) == "ATOM" or line.substr(0, 6) == "HETATM") {
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
  confFile << "size_x = " << sizex + 10 << std::endl;
  confFile << "size_y = " << sizey + 10 << std::endl;
  confFile << "size_z = " << sizez + 10 << std::endl;
  confFile.close();
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
  std::string gromacsPath = reader.Get("paths", "gromacs", "gmx");
  std::string mdpPath = reader.Get("paths", "mdp", "");
  std::string forcefieldPath = reader.Get("paths", "forcefieldpath", "");
  std::string forcefield = reader.Get("GROMACS", "forcefield", "");
  std::string water = reader.Get("GROMACS", "water", "");
  std::string boundingboxtype = reader.Get("GROMACS", "bt", "");
  float boxsize = reader.GetReal("GROMACS", "bt", 1.0);
  float clustercutoff = reader.GetReal("GROMACS", "clustercutoff", 0.12);

  std::cout << gromacsPath << std::endl;
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
  // Generate config files
  for (auto s : receptors) {
    prepareConfig(s);
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
