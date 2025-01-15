#ifndef IO_H
#define IO_H

#include <string>
#include <vector>
#include <fstream>
#include "particles.hpp"

configuration ReadTrimCFG(std::string input);

void WriteTrimCFG(const configuration& cfg, std::string output);

std::ofstream MakeObsFile(std::vector <std::string>& observables, std::string output);

void WriteObs(const configuration& cfg, const configuration& cfg0, 
              int t, int cycle, std::vector <std::string>& observables, 
              std::ofstream& log_obs);

std::vector <std::pair <int,int>> GetLogspacedSnapshots(int cycles, int tau, int tw, int n_log);

std::vector <int> GetLinspacedSnapshots(int cycles, int tau, int tw, int n_lin);

#endif