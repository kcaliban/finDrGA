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


float VinaInstance::calculateBindingAffinity(int exhaustiveness,
                                              int energy_range) {
  // Create command
  std::string command;
  command.append(vinaPath);
  command.append(" --config ");
  command.append(receptor);
  command.append("_conf");
  command.append(" --exhaustiveness ");
  command.append(std::to_string(exhaustiveness));
  command.append(" --receptor ");
  command.append(receptor);
  // Receptor is already in qt form (in the case of prepared input)
  if (receptor.substr(receptor.size() - 2, 2) != "qt") {
    command.append("qt");
  }
  command.append(" --ligand ");
  command.append(ligand);
  command.append("qt");
  command.append(" --energy_range ");
  command.append(std::to_string(energy_range));

  std::string out = "Docking against " + receptor + "...";
  debugPrint(out.c_str());

  // Execute command and stream using popen
  std::string vinaOutput;
  FILE * vinaOutptStream = popen(command.c_str(), "r");
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
