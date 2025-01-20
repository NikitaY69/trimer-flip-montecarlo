#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include "particles.hpp"

/**
 * @brief Runs the Monte Carlo simulation.
 * 
 * @param cfg Initial configuration.
 * @param T Temperature.
 * @param tau Number of Monte Carlo steps.
 * @param cycles Number of cycles.
 * @param tw Waiting time.
 * @param p_flip Probability of flipping.
 * @param observables List of observables to compute.
 * @param out Output directory.
 * @param n_log Number of log-spaced points.
 * @param n_lin Number of linear-spaced points.
 */
void MonteCarloRun(configuration& cfg, double T, int tau, int cycles, int tw, double p_flip, 
        std::vector <std::string>& observables, std::string& out, int n_log, int n_lin);

/**
 * @brief Tries displacing one particle.
 * 
 * @param cfg Current configuration.
 * @param j Index of the particle to displace.
 * @param T Temperature.
 */
void TryDisp(configuration& cfg, int j, double T);

/**
 * @brief Tries swapping two particles' diameters.
 * 
 * @param cfg Current configuration.
 * @param j Index of the particle to swap.
 * @param T Temperature.
 */
void TryFlip(configuration& cfg, int j, double T);

/**
 * @brief Computes observables without running the simulation.
 * 
 * @param tau Number of Monte Carlo steps.
 * @param cycles Number of cycles.
 * @param tw Waiting time.
 * @param observables List of observables to compute.
 * @param out Output directory.
 * @param n_log Number of log-spaced points.
 */
void ComputeObservables(int tau, int cycles, int tw, 
        std::vector <std::string>& observables, std::string& out, int n_log);

#endif // SIMULATION_H