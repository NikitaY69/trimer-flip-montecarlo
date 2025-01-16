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
    std::string rootdir;
    std::vector <std::string> observables;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Produce help message")
        ("init", po::value<std::string>(&input), "Path to initial configuration file")
        ("params", po::value<std::string>()->required(), "Path to JSON file for simulation parameters")
        ("observables", po::value<std::vector<std::string>>(&observables)->multitoken(),
                        "List of observables to compute (e.g., MSD Fs U; separated by spaces)");

    // Parse the command-line arguments
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        // Handle the help option
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        po::notify(vm);
    } catch (const po::error &ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    // Loading params from json file
    json params;
    std::ifstream input_file(vm["params"].as<std::string>());

    // Check if the file is open
    if (!input_file.is_open()) {
        std::cerr << "Error opening file.\n";
        return 1;
    }

    // Parse the JSON file
    input_file >> params;
    rootdir = params["rootdir"];
    N = params["N"];
    T = params["T"];
    tau = params["tau"];
    tw = params["tw"];
    cycles = params["cycles"];
    logPoints = params["logPoints"];
    linPoints = params["linPoints"];
    p_flip = params["p_flip"];

    // Recalculating size
    Size = pow(N/density, 1./3.);

    // Creating outdir if not existing
    fs::path out_path = rootdir;
    if(!fs::is_directory(out_path)){
        fs::create_directory(rootdir);
    }

    std::ofstream file(rootdir + "params.json");
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
    MC(initconf, T, tau, cycles, tw, p_flip, observables, rootdir, logPoints, linPoints); 
    double time_elapsed = time(NULL) - t0;
    std::cout << "================================================================" 
              << std::endl;
    std::cout << "Time taken: " << std::endl;
    std::cout << time_elapsed << " seconds" << std::endl;
    std::cout << time_elapsed/60 << " minutes" << std::endl;
    std::cout << time_elapsed/3600 << " hours" << std::endl;

    return 0;
}