#include "VinaInstance.h"

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

void VinaInstance::generatePDBQT() {
  // Generate receptor
  char cmd[3000];
  strcpy(cmd, pythonShPath.c_str());
  strcpy(cmd, " ");
  strcat(cmd, mgltoolstilitiesPath.c_str());
  strcat(cmd, "/prepare_receptor4.py -r ");
  strcat(cmd, receptor.c_str());
  strcat(cmd, " -A bonds_hydrogens -U nphs -o ");
  strcat(cmd, receptor.c_str());
  strcat(cmd, "qt");
  strcat(cmd, " >/dev/null"); // Hides console output
  int success = system(cmd);

  if (success == -1) {
    std::cout << "ERROR: Could not generate pdbqt file from receptor" << std::endl;
    exit(-1);
  }

  memset(cmd, 0, sizeof(cmd));

  // Generate ligand
  strcpy(cmd, pythonShPath.c_str());
  strcpy(cmd, " ");
  strcat(cmd, mgltoolstilitiesPath.c_str());
  strcat(cmd, "/prepare_ligand4.py -l ");
  strcat(cmd, ligand.c_str());
  strcat(cmd, " -A bonds_hydrogens -U nphs -o ");
  strcat(cmd, ligand.c_str());
  strcat(cmd, "qt");
  strcat(cmd, " >/dev/null");
  success = system(cmd);

  if (success == -1) {
    std::cout << "ERROR: Could not generate pdbqt file from receptor" << std::endl;
    exit(-1);
  }
}

float VinaInstance::calculateBindingAffinity(int exhaustiveness,
                                              int energy_range) {
  // Create command
  char cmd[3000];
  strcpy(cmd, vinaPath.c_str());
  strcat(cmd, "/vina");
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

  std::cout << "\t\t\tDocking of: " << ligand << std::endl;

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

  /* For debugging
  std::cout << vinaOutput << std::endl;

  std::cout << affinityMatch.size() << std::endl;
  std::cout << affinityMatch.str(0) << std::endl;
  std::cout << "Match: " << affinityMatch.str(1) << ", " << stof(affinityMatch.str(1)) << std::endl;
  */
  return stof(affinityMatch.str(1));
}
