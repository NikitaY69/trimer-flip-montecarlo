/**
 * @file utils.hpp
 * @brief Utility functions for parsing input, reading/writing configuration files,
 *        handling simulation parameters, and managing observables.
 *
 * This module provides helper functions to parse command-line arguments,
 * read/write JSON configuration files, handle trimer configurations, and
 * manage simulation observables. It is designed to streamline the setup and
 * execution of simulations involving particles and their properties.
 */

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <fstream>
#include "particles.hpp"

/**
 * @brief Parses command line arguments.
 * 
 * @param argc Number of arguments.
 * @param argv Array of argument strings.
 * @param input Path to initial configuration file.
 * @param params Path to JSON file for simulation parameters.
 * @param observables List of observables to compute.
 * @param seed Random seed.
 * @return true if parsing was successful, false otherwise.
 */
bool ParseCMDLine(int argc, const char* argv[],
                        std::string& input,
                        std::string& params,
                        std::vector<std::string>& observables,
                        int& seed);

/**
 * @brief Reads parameters from a JSON file.
 * 
 * @param params_path Path to the JSON file.
 * @param rootdir Root directory for output files.
 * @param N Number of particles.
 * @param T Temperature.
 * @param tau Number of Monte Carlo steps.
 * @param tw Waiting time.
 * @param cycles Number of cycles.
 * @param logPoints Number of log-spaced points.
 * @param linPoints Number of linear-spaced points.
 * @param p_flip Probability of flipping.
 * @return true if reading was successful, false otherwise.
 */
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

/**
 * @brief Reads trimer configuration from a file.
 * 
 * @param input Path to the input file.
 * @return The configuration read from the file.
 */
configuration ReadTrimCFG(std::string input);

/**
 * @brief Writes trimer configuration to a file.
 * 
 * @param cfg The configuration to write.
 * @param output Path to the output file.
 */
void WriteTrimCFG(const configuration& cfg, std::string output);

/**
 * @brief Creates output directory for storing results.
 * 
 * @param rootdir Root directory for output files.
 * @param params_path Path to the parameter file (used for naming).
 */
void MakeOutDir(std::string rootdir, std::string params_path);

/**
 * @brief Creates an output file for storing observable data.
 * 
 * @param observables List of observables to log.
 * @param output Path to the output file.
 * @return An open output file stream for logging observables.
 */
std::ofstream MakeObsFile(std::vector <std::string>& observables, std::string output);

/**
 * @brief Writes observables at a specific timestep.
 * 
 * @param cfg Current configuration.
 * @param cfg0 Initial configuration.
 * @param t Current timestep.
 * @param cycle Current cycle.
 * @param observables List of observables.
 * @param log_obs Output file stream for logging observables.
 */
void WriteObs(const configuration& cfg, const configuration& cfg0, 
              int t, int cycle, std::vector <std::string>& observables, 
              std::ofstream& log_obs);

/**
 * @brief Gets log-spaced snapshots.
 * 
 * @param cycles Number of cycles.
 * @param tau Number of Monte Carlo steps.
 * @param tw Waiting time.
 * @param n_log Number of log-spaced points.
 * @return Vector of pairs of log-spaced snapshots.
 */
std::vector <std::pair <int,int>> GetLogspacedSnapshots(int cycles, int tau, int tw, int n_log);

/**
 * @brief Gets linear-spaced snapshots.
 * 
 * @param cycles Number of cycles.
 * @param tau Number of Monte Carlo steps.
 * @param tw Waiting time.
 * @param n_lin Number of linear-spaced points.
 * @return Vector of linear-spaced snapshots.
 */
std::vector <int> GetLinspacedSnapshots(int cycles, int tau, int tw, int n_lin);

#endif // UTILS_H
