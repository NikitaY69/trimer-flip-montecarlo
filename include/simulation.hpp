#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include "particles.hpp"

void MonteCarloRun(configuration& cfg, double T, int tau, int cycles, int tw, double p_flip, 
        std::vector <std::string>& observables, std::string& out, int n_log, int n_lin);

void TryDisp(configuration& cfg, int j, double T);

void TryFlip(configuration& cfg, int j, double T);

#endif