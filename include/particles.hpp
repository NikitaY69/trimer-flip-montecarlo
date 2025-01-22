/**
 * @file particles.hpp
 * @brief Structures and functions for particle configurations and interactions.
 *
 * This module defines the `configuration` structure to represent the state of particles
 * in a simulation, including their positions, bonded interactions, and Verlet lists.
 * It also provides utility functions for handling periodic boundary conditions and
 * managing particle coordinates.
 */

#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include "globals.hpp"

/**
 * @brief Structure to keep track of the evolution of configurations.
 */
struct configuration {
    std::vector<double> X;      ///< Particles' X coordinates inside main box
    std::vector<double> Y;      ///< Particles' Y coordinates inside main box
    std::vector<double> Z;      ///< Particles' Z coordinates inside main box
    std::vector<double> Xfull;  ///< Particles' real X coordinates (for dynamical purposes)
    std::vector<double> Yfull;  ///< Particles' real Y coordinates (for dynamical purposes)
    std::vector<double> Zfull;  ///< Particles' real Z coordinates (for dynamical purposes)
    std::vector<double> X0;     ///< Particles' X coordinates at last neighbors update
    std::vector<double> Y0;     ///< Particles' Y coordinates at last neighbors update
    std::vector<double> Z0;     ///< Particles' Z coordinates at last neighbors update
    std::vector<int> S;         ///< Particles' types (NOT DIAMETERS)
    std::vector<std::vector<int>> neighbours_list; ///< Verlet lists for each particle
    std::vector<std::vector<int>> bonded_neighbours; ///< Bonded particles for each particle
    double XCM; ///< Center of mass X coordinate
    double YCM; ///< Center of mass Y coordinate
    double ZCM; ///< Center of mass Z coordinate
    
    /**
     * @brief Constructor to initialize vectors.
     */
    configuration() 
        : X(N), Y(N), Z(N), Xfull(N), Yfull(N), Zfull(N), 
          X0(N), Y0(N), Z0(N), S(N), neighbours_list(N), bonded_neighbours(N), 
          XCM(0), YCM(0), ZCM(0) {}

    /**
     * @brief Method to calculate center of mass coordinates.
     */
    void UpdateCM_coord();

    /**
     * @brief Method to update the verlet lists for each particle.
     */
    void UpdateNL();

    /**
     * @brief Method to retrieve bonded particles for all particles.
     */
    void GetBonds();

    /**
     * @brief Method to check whether or not to update the neighbors list.
     */
    void CheckNL();
};

/**
 * @brief Calculates difference of a and b while applying periodic boundary conditions.
 * 
 * @param a First coordinate
 * @param b Second coordinate
 * @return Adjusted coordinate difference
 */
double bcs(double a, double b);

/**
 * @brief Shifts coordinate inside main box.
 * 
 * @param a Coordinate to be shifted
 * @return Shifted coordinate
 */
double Pshift(double a);

#endif // PARTICLES_H