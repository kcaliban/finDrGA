#include "VinaInstance.h"

std::string VinaInstance::logStr() {
  if (log) {
    return " >> " + ligand.substr(0,
        ligand.size() - 15) + "/VINAINSTLOG" + " 2>&1";
  }
  return " > /dev/null 2>&1";
}

void VinaInstance::debugPrint(const char * str) {
  if (debug) {
    std::cout << "\033[1;32mINFO (DOCKING, " << ligand << "): " << str << "\033[0m" << std::endl;
  }
}

void VinaInstance::errorPrint(const char * str) {
  std::cout << "\033[1;31mERROR (DOCKING): " << str << std::endl;
  std::cout << "Ligand: " << ligand << "\033[0m" << std::endl;

  exit(-1);
}

void VinaInstance::generateConf() {
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

  std::string outfile = workDir;
  /*
  std::cout << "Workdir @ generateConf" << std::endl;
  std::cout << workDir << std::endl;
  std::cout << "---------------------" << std::endl;
  */
  outfile.append("/");
  outfile.append("conf");

  std::ofstream confFile;
  confFile.open(outfile.c_str(), std::ios::trunc);
  if (!confFile) {
    errorPrint("Could not open config file!");
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

void VinaInstance::generatePDBQT() {
  debugPrint("Generating pdbqts...");
  // Generate receptor PDBQT
  char cmd[3000];
  strcpy(cmd, pythonShPath.c_str());
  strcat(cmd, " ");
  strcat(cmd, mgltoolstilitiesPath.c_str());
  strcat(cmd, "/prepare_receptor4.py -r ");
  strcat(cmd, receptor.c_str());
  strcat(cmd, " -A bonds_hydrogens -U nphs -o ");
  strcat(cmd, receptor.c_str());
  strcat(cmd, "qt");
  strcat(cmd, logStr().c_str());
  int success = system(cmd);

  if (success != 0) {
    errorPrint("Could not generate pdbqt file from receptor");
  }

  memset(cmd, 0, sizeof(cmd));

  // Generate ligand PDBQT
  strcpy(cmd, pythonShPath.c_str());
  strcat(cmd, " ");
  strcat(cmd, mgltoolstilitiesPath.c_str());
  strcat(cmd, "/prepare_ligand4.py -l ");
  strcat(cmd, ligand.c_str());
  strcat(cmd, " -A bonds_hydrogens -U nphs -o ");
  strcat(cmd, ligand.c_str());
  strcat(cmd, "qt");
  strcat(cmd, logStr().c_str());
  success = system(cmd);

  if (success != 0) {
    errorPrint("Could not generate pdbqt file from receptor");
  }
}

float VinaInstance::calculateBindingAffinity(int exhaustiveness,
                                              int energy_range) {
  // Create command
  char cmd[3000];
  strcpy(cmd, vinaPath.c_str());
  strcat(cmd, " --config ");
  strcat(cmd, workDir.c_str());
  strcat(cmd, "/conf");
  strcat(cmd, " --exhaustiveness ");
  strcat(cmd, (std::to_string(exhaustiveness)).c_str());
  strcat(cmd, " --receptor ");
  strcat(cmd, receptor.c_str());
  strcat(cmd, "qt");
  strcat(cmd, " --ligand ");
  strcat(cmd, ligand.c_str());
  strcat(cmd, "qt");
  strcat(cmd, " --energy_range ");
  strcat(cmd, (std::to_string(energy_range)).c_str());

  debugPrint("Docking...");

  // Execute command and stream using popen
  std::string vinaOutput;
  FILE * vinaOutptStream = popen(cmd, "r");
  char buff[1024];
  while (fgets(buff, 1024, vinaOutptStream)) {
    vinaOutput += buff;
  }
  pclose(vinaOutptStream);

  // Use a regex to find the best affinity
  std::regex affinityRegEx("\n   1[ ]*([-.0-9]+)");
  std::smatch affinityMatch;
  std::regex_search(vinaOutput, affinityMatch, affinityRegEx);

  return stof(affinityMatch.str(1));
}
