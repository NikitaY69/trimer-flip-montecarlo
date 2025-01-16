#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <ctime>
#include <experimental/filesystem>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include "globals.hpp"
#include "particles.hpp"
#include "utils.hpp"
#include "montecarlo.hpp"

namespace fs = std::experimental::filesystem;
namespace po = boost::program_options;
using json = nlohmann::json;

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
    srand(time(NULL)*1.0);

    // Define the command-line options
    std::string input;
    std::string outdir;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input", po::value<std::string>(&input)->required(), "set input file")
        ("outdir", po::value<std::string>(&outdir)->required(), "set out directory")
        ("N", po::value<int>(&N)->default_value(N), "set system size")
        ("T", po::value<double>(&T)->default_value(T), "set temperature")
        ("tau", po::value<int>(&tau)->default_value(tau), "set single-run time")
        ("tw", po::value<int>(&tw)->default_value(tw), "set waiting time")
        ("cycles", po::value<int>(&cycles)->default_value(cycles), "set number of cycles")
        ("lin", po::value<int>(&linPoints)->default_value(linPoints), "set number of lin-spaced snapshots")
        ("log", po::value<int>(&logPoints)->default_value(logPoints), "set number of log-spaced snapshots")
        ("p_flip", po::value<double>(&p_flip)->default_value(p_flip), "set flip-attempt probability")
        ("MSD", "Flag to compute MSD")
        ("Fs", "Flag to compute Fs")
        ("U", "Flag to compute U");

    // Parse the command-line arguments
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const po::error &ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    // Handle the help option
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    // Parsing the observables in order of appearance
    std::vector < std::string > observables;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--MSD") {
            observables.push_back("MSD");
        } else if (arg == "--Fs") {
            observables.push_back("Fs");
        } else if (arg == "--U") {
            observables.push_back("U");
        }
    }

    // Recalculating user-defined parameters
    Size = pow(N/density, 1./3.);

    // Creating outdir if not existing
    fs::path out_path = outdir;
    if(!fs::is_directory(out_path)){
        fs::create_directory(outdir);
    }
    
     // Writing params.json file
    json params;
    params["rootdir"] = outdir;
    params["N"] = N;
    params["T"] = T;
    params["tau"] = tau;
    params["tw"] = tw;
    params["cycles"] = cycles;
    params["logPoints"] = logPoints;
    params["linPoints"] = linPoints;
    params["p_flip"] = p_flip;

    std::ofstream file(outdir + "params.json");
    if (file.is_open()) {
        file << params.dump(4); // Pretty-print JSON with 4 spaces of indentation
        file.close();
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }

    // Read init config
    configuration initconf;
    initconf = ReadTrimCFG(input);

    // Do simulation with timer
    double t0 = time(NULL); // Timer
    MC(initconf, T, tau, cycles, tw, p_flip, 
       observables, outdir, logPoints, linPoints); 
    std::cout << "Time taken: " << (time(NULL) - t0) << "s" << std::endl; 
    std::cout << "Done" << std::endl;

    return 0;
}