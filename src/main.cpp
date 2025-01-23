#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include "globals.hpp"
#include "particles.hpp"
#include "utils.hpp"
#include "simulation.hpp"

// Default run parameters
int N = 5;
double Size = pow(N/density, 1/3.);
double T = 2.0; 
int tau = 100000;
int tw = 1;
int cycles = 1;
int linPoints = 50;
int logPoints = 50;
double p_flip = 0.2;

//-----------------------------------------------------------------------------
//  main.cpp
int main(int argc, const char * argv[]) {

    // Define the command-line options
    std::string input;
    std::string params_path;
    std::string rootdir;
    std::vector <std::string> observables;
    bool norun;
    
    // Parse command line arguments
    if (!ParseCMDLine(argc, argv, input, params_path, observables)){
        return 1;
    };

    // Loading params from json file
    if (!ReadJSONParams(params_path, rootdir, N, T, tau, tw, cycles, logPoints, linPoints, p_flip)){
        return 1;
    }

    // Random number seed
    srand(time(NULL));

    // Recalculating size
    Size = pow(N/density, 1./3.);

    // if (fs::exists(target_path)) {
    //     // pass
    // } else {
    // Copy the file to the target directory
    
    // Setting run mode
    (input == "") ? norun = true : norun = false;

    double t0 = time(NULL); // Timer

    if (norun){
        // Compute observables
        ComputeObservables(tau, cycles, tw, observables, rootdir, logPoints);
    } else{
        // Read init config
        configuration initconf;
        initconf = ReadTrimCFG(input);
        // Make outdir and copy json file
        MakeOutDir(rootdir, params_path);
        // Do simulation
        MonteCarloRun(initconf, T, tau, cycles, tw, p_flip, observables, rootdir, logPoints, linPoints, true); 
    }
    
    double time_elapsed = time(NULL) - t0;
    std::cout << "Time taken: " << std::endl;
    std::cout << time_elapsed << " seconds" << std::endl;
    std::cout << time_elapsed/60 << " minutes" << std::endl;
    std::cout << time_elapsed/3600 << " hours" << std::endl;

    return 0;
}