/**
 * @file globals.hpp
 * @brief Global variables and macros used across the simulation.
 *
 * This module defines shared global variables and constants, such as the number of particles,
 * simulation box size, and particle diameters. It also includes a macro for generating random
 * numbers within the range [0, 1].
 *
 * The constants and globals defined here are used to configure and manage simulation parameters
 * at a global scope.
 */

#ifndef GLOBALS_H
#define GLOBALS_H

// Libraries
#include <cstdlib> // For dynamic memory and random numbers

// Shared variables

/**
 * @brief Number of particles.
 */
extern int N;

/**
 * @brief Size of the simulation box.
 */
extern double Size;

/**
 * @brief Maximum diameter of particles.
 */
const double sigmaMax = 1.1;

/**
 * @brief Array of particle diameters.
 */
const double diameters[3] = {0.9, 1.0, 1.1};

/**
 * @brief Generates a random number between 0 and 1.
 */
#define ranf() \
    ((double)rand()/(1.0+RAND_MAX))

#endif // GLOBALS_H