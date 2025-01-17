#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <fstream>
#include "particles.hpp"

bool ParseCMDLine(int argc, const char* argv[],
                        std::string& input,
                        std::string& params,
                        std::vector<std::string>& observables,
                        int& seed);

bool ReadJSONParams(const std::string& params_path, 
                    std::string& rootdir,
                    int& N, 
                    double& T,
                    int& tau,
                    int& tw,
                    int& cycles,
                    int& logPoints,
                    int& linPoints,
                    double& p_flip);
                    
configuration ReadTrimCFG(std::string input);

void WriteTrimCFG(const configuration& cfg, std::string output);

std::ofstream MakeObsFile(std::vector <std::string>& observables, std::string output);

void WriteObs(const configuration& cfg, const configuration& cfg0, 
              int t, int cycle, std::vector <std::string>& observables, 
              std::ofstream& log_obs);

std::vector <std::pair <int,int>> GetLogspacedSnapshots(int cycles, int tau, int tw, int n_log);

std::vector <int> GetLinspacedSnapshots(int cycles, int tau, int tw, int n_lin);

#endif