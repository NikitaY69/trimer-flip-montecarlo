/**
 * @file simulation.hpp
 * @brief Core simulation functions for Monte Carlo methods.
 *
 * This module provides the main functions for running Monte Carlo simulations,
 * including particle displacement, diameter swapping, and observable computations.
 * These functions are the building blocks for simulations involving particle
 * interactions and dynamics at various temperatures and configurations.
 */

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
 * @param progress_bar True/False statement to output progress_bar
 */
void MonteCarloRun(configuration& cfg, double T, int tau, int cycles, int tw, double p_flip, 
        std::vector <std::string>& observables, std::string& out, int n_log, int n_lin, bool progress_bar);

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