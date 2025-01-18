#ifndef GLOBALS_H
#define GLOBALS_H

// Libraries
#include <cstdlib> // For dynamic memory and random numbers

// Shared variables
const double density = 1.2; // Density of the system
extern int N; // Number of particles
extern double Size; // Size of simulation box
const double sigmaMax = 1.1; // Maximum diameter of particles
const double diameters[3] = {0.9, 1.0, 1.1};

//  Random number between 0 and 1
#define ranf() \
    ((double)rand()/(1.0+RAND_MAX))
    
#endif