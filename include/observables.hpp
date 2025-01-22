/**
 * @file obs.hpp
 * @brief Functions for calculating observables in particle simulations.
 *
 * This module provides functions to calculate various observables and potentials
 * used in particle-based simulations. These include pairwise potentials like WCA
 * and FENE, as well as system-level metrics such as total energy, mean square
 * displacement, and intermediate self-scattering functions.
 */

#ifndef OBS_H
#define OBS_H

#include "particles.hpp"

/**
 * @brief Calculates the pairwise WCA potential between two particles.
 * 
 * @param x1 X coordinate of the first particle.
 * @param y1 Y coordinate of the first particle.
 * @param z1 Z coordinate of the first particle.
 * @param s1 Diameter of the first particle.
 * @param x2 X coordinate of the second particle.
 * @param y2 Y coordinate of the second particle.
 * @param z2 Z coordinate of the second particle.
 * @param s2 Diameter of the second particle.
 * @return The WCA potential between the two particles.
 */
double WCAPair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2);

/**
 * @brief Calculates the pairwise FENE potential between two particles.
 * 
 * @param x1 X coordinate of the first particle.
 * @param y1 Y coordinate of the first particle.
 * @param z1 Z coordinate of the first particle.
 * @param s1 Diameter of the first particle.
 * @param x2 X coordinate of the second particle.
 * @param y2 Y coordinate of the second particle.
 * @param z2 Z coordinate of the second particle.
 * @param s2 Diameter of the second particle.
 * @return The FENE potential between the two particles.
 */
double FENEPair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2);

/**
 * @brief Calculates the potential associated with a particle.
 * 
 * @param cfg Current configuration.
 * @param j Index of the particle.
 * @return The potential associated with the particle.
 */
double V(const configuration& cfg, int j);

/**
 * @brief Calculates the total system energy.
 * 
 * @param cfg Current configuration.
 * @return The total system energy.
 */
double VTotal(const configuration& cfg);

/**
 * @brief Calculates the average mean square displacement.
 * 
 * @param cfg Current configuration.
 * @param cfg0 Initial configuration.
 * @return The average mean square displacement.
 */
double MSD(const configuration& cfg, const configuration& cfg0);

/**
 * @brief Calculates the intermediate self-scattering function.
 * 
 * @param cfg Current configuration.
 * @param cfg0 Initial configuration.
 * @return The intermediate self-scattering function.
 */
double FS(const configuration& cfg, const configuration& cfg0);

#endif // OBS_H