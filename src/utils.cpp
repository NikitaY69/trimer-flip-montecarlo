#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/json.hpp>
#include <boost/filesystem.hpp>
#include "globals.hpp"
#include "utils.hpp"
#include "observables.hpp"

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace json = boost::json;

// Parse command line arguments
bool ParseCMDLine(int argc, const char* argv[],
                        std::string& input,
                        std::string& params,
                        std::vector<std::string>& observables,
                        int& seed) {
    // Define the command-line options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Produce help message")
        ("init", po::value<std::string>(&input)->default_value(""), "Path to initial configuration file")
        ("params", po::value<std::string>(&params)->required(), "Path to JSON file for simulation parameters")
        ("observables", po::value<std::vector<std::string>>(&observables)->multitoken(),
                        "List of observables to compute (e.g., MSD Fs U; separated by spaces)")
        ("seed", po::value<int>(&seed)->default_value(static_cast<int>(time(NULL))), "random seed");

    // Parse the command-line arguments
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        // Handle the help option
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return false; // Indicate that help was requested
        }

        po::notify(vm);
    } catch (const po::error& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return false;
    }

    return true; // Successfully parsed arguments
}

// Read params from JSON file
bool ReadJSONParams(const std::string& params_path, 
                    std::string& rootdir,
                    int& N, 
                    double& T,
                    int& tau,
                    int& tw,
                    int& cycles,
                    int& logPoints,
                    int& linPoints,
                    double& p_flip) {
    std::ifstream json_file(params_path);

    // Check if the file is open
    if (!json_file.is_open()) {
        std::cerr << "Error opening JSON file.\n";
        return false;
    }

    // Read and parse the JSON file
    std::stringstream buffer;
    buffer << json_file.rdbuf();
    json::value parsed_json = json::parse(buffer.str());

    // Access JSON fields
    const auto& obj = parsed_json.as_object();
    rootdir = obj["rootdir"].as_string().c_str();
    N = obj["N"].as_int64();
    T = obj["T"].as_double();
    tau = obj["tau"].as_int64();
    tw = obj["tw"].as_int64();
    cycles = obj["cycles"].as_int64();
    logPoints = obj["logPoints"].as_int64();
    linPoints = obj["linPoints"].as_int64();
    p_flip = obj["p_flip"].as_double();

    return true;
}

// Read trimer configs
configuration ReadTrimCFG(std::string input){
    configuration C;
    int type;
    std::string line;
    std::ifstream input_file(input);
    if (!input_file.is_open()){
        throw std::runtime_error("Could not open file: " + input);
    }
    int i = 0; // particle index
    std::vector<std::vector<double>> cfg; // array of configurations
    while (std::getline(input_file, line)){
        double value;
        std::stringstream ss(line);

        cfg.push_back(std::vector<double>());
        while (ss >> value){
            cfg[i].push_back(value);
        }
        // mol_index[i] = cfg[i][0];
        // type = cfg[i][1]-1; 
        // C.S[i] = (diameters[type]);
        if (cfg[i].size() == 5){
            C.S[i] = cfg[i][1];
            C.Xfull[i] = cfg[i][2]; C.Yfull[i] = cfg[i][3]; C.Zfull[i] = cfg[i][4];
        } else{
            C.S[i] = cfg[i][0];
            C.Xfull[i] = cfg[i][1]; C.Yfull[i] = cfg[i][2]; C.Zfull[i] = cfg[i][3];
        }
        
        C.X[i] = Pshift(C.Xfull[i]); C.X0[i] = C.X[i]; 
        C.Y[i] = Pshift(C.Yfull[i]); C.Y0[i] = C.Y[i];
        C.Z[i] = Pshift(C.Zfull[i]); C.Z0[i] = C.Z[i];
        i++;}
    input_file.close();
    return C;
}

// Write trimer configs
void WriteTrimCFG(const configuration& cfg, std::string output){
    std::ofstream log_cfg;
    log_cfg.open(output);
    log_cfg << std::scientific << std::setprecision(8);
    for (int i = 0; i<N; i++){
        log_cfg << cfg.S[i] << " " << cfg.Xfull[i] << " " << cfg.Yfull[i] << " " << cfg.Zfull[i] << std::endl;
    }
    log_cfg.close();
}

// Make output directory
void MakeOutDir(std::string rootdir, std::string params_path) {
    fs::path rootdir_path = rootdir;
    if (!fs::is_directory(rootdir_path)) {
        fs::create_directory(rootdir);
    }

    // Configs directory
    std::string out_cfg = rootdir + "configs/";
    fs::create_directory(out_cfg);

    // Copy the params file into the root directory
    fs::path json_file(params_path);
    fs::path target_path = rootdir_path / json_file.filename();
    try {
        fs::copy_file(json_file, target_path, fs::copy_options::overwrite_existing);
    } catch (const fs::filesystem_error& e) {
        std::cerr << e.what() << std::endl;
    }
}

// Create observables file
std::ofstream MakeObsFile(std::vector <std::string>& observables, std::string output){
    std::ofstream log_obs; 
    log_obs.open(output);
    log_obs << "t" << " " << "cycle";
    for (const std::string obs: observables){
        log_obs << " " << obs;
    } log_obs << std::endl;
    log_obs << std::scientific << std::setprecision(8);
    return log_obs;
}

// Write observables at specific timestep
void WriteObs(const configuration& cfg, const configuration& cfg0, 
              int t, int cycle, std::vector <std::string>& observables, 
              std::ofstream& log_obs){
    
    log_obs << t << " " << cycle;
    for (std::string obs: observables){
        log_obs << " ";
        (obs == "U")   ? log_obs << VTotal(cfg)/(2*N) : 
        (obs == "MSD") ? log_obs << MSD(cfg, cfg0) : 
                         log_obs << FS(cfg, cfg0);
    } log_obs << std::endl;
}

std::vector <std::pair <int,int>> GetLogspacedSnapshots(int cycles, int tau, int tw, int n_log){
    std::vector < std::pair <int, int>> pairs;
    std::vector <int> samplePoints, twPoints;
    double exponents = log10(tau)/(n_log-1);

    for(int c=0; c<cycles; c++){
        for (int x = 0; x < n_log; x++){
            int value = tw*c + floor(pow(10,exponents*(x)));
            std::pair <int,int> p = {value, c};
            int f = std::count(pairs.begin(), pairs.end(), p);
            if(f==0){
                pairs.emplace_back(value, c);
            // this if condition is actually relevent because of the floor function
            }
        }
    }

    // Sorting
    std::sort(pairs.begin(), pairs.end());
    return pairs;

}

std::vector <int> GetLinspacedSnapshots(int cycles, int tau, int tw, int n_lin){
    std::vector <int> linpoints;
    for (int c=0; c<cycles; c++){
        for (int k=1; k<=n_lin; k++){
            linpoints.push_back(tw*c+(tau/(n_lin))*k);
        }
    } std::sort(linpoints.begin(), linpoints.end());
    return linpoints;
}