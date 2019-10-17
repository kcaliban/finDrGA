/* Copyright 2019 iGEM Team Freiburg 2019 */
#include "VinaInstance.h"

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
  command.append(" --cpu 1");
  command.append(" --energy_range ");
  command.append(std::to_string(energy_range));
  command.append(" --out ");
  std::string outName = receptor.substr(receptor.find_last_of("/") + 1,
                                        receptor.size() -
                                        receptor.find_last_of("/") - 1);
  command.append(ligand + outName);

  command.append(" --log ");
  command.append(ligand + outName + "VINALOG");

  info->infoMsg("(VINA) Docking " + ligand + " against: " + receptor);

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
  try {
    std::regex_search(vinaOutput, affinityMatch, affinityRegEx);
  } catch (std::regex_error &e) {
    info->errorMsg("Error in VINA output. Output:\n" + vinaOutput, true);
  }

  return stof(affinityMatch.str(1));
}
