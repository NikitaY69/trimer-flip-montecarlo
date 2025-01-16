#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include <experimental/filesystem>
#include "globals.hpp"
#include "particles.hpp"
#include "utils.hpp"
#include "simulation.hpp"

namespace fs = std::experimental::filesystem;

// Default run parameters
const double density = 1.2;
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

    // Random number seed
    srand(40);

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

    // Recalculating size
    Size = pow(N/density, 1./3.);

    // Creating outdir if not existing
    fs::path rootdir_path = rootdir;
    if(!fs::is_directory(rootdir_path)){
        fs::create_directory(rootdir);
    }

    // Checking if the json file is present in rootdir
    fs::path json_file(params_path);
    fs::path target_path = rootdir_path / json_file.filename();

    // if (fs::exists(target_path)) {
    //     // pass
    // } else {
    // Copy the file to the target directory
    try {
        fs::copy(json_file, target_path, fs::copy_options::overwrite_existing);
    } catch (const fs::filesystem_error& e) {
        std::cerr << e.what() << std::endl;
    }
    
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

        // Do simulation
        MonteCarloRun(initconf, T, tau, cycles, tw, p_flip, observables, rootdir, logPoints, linPoints); 
    }
    
    double time_elapsed = time(NULL) - t0;
    std::cout << "Time taken: " << std::endl;
    std::cout << time_elapsed << " seconds" << std::endl;
    std::cout << time_elapsed/60 << " minutes" << std::endl;
    std::cout << time_elapsed/3600 << " hours" << std::endl;

    return 0;
}