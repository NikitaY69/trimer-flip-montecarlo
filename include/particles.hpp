#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include "globals.hpp"

// Structure to keep track of the evolution of configurations
struct configuration {
    std::vector <double> X, Y, Z, Xfull, Yfull, Zfull, X0, Y0, Z0;
    std::vector <int> S;
    // X,Y,Z: size N vectors containing particles coordinates inside main box
    // Xfull, Yfull, Zfull: size N vectors containing particles 
    // X0, Y0, Z0: size N vectors containing particles coordinates at last neighbours update
    // S: size N vector containing particles diameters
    std::vector < std::vector<int> > NL, BN; // verlet lists and bonded particles for each particle
    double XCM, YCM, ZCM; // center of mass coordinates
    
    // Constructor to initialize vectors
    configuration(): X(N), Y(N), Z(N), Xfull(N), Yfull(N), Zfull(N), 
                     X0(N), Y0(N), Z0(N), S(N), // 
                     NL(N), BN(N),
                     XCM(0), YCM(0), ZCM(0){};
    // Method to calculate center of mass coordinates
    void UpdateCM_coord();
    // Method to update the verlet lists for each particle
    void UpdateNL(), GetBonds();
    // Method to check whether or not to update the neighbours list
    void CheckNL();
};

double bcs(double a, double b), Pshift(double a);

#endif