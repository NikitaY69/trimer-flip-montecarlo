#include <catch2/catch.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <boost/filesystem.hpp>
#include "globals.hpp"
#include "utils.hpp"
#include "simulation.hpp"

namespace fs = boost::filesystem;

// Function to check if two files are identical
bool AreFilesIdentical(const std::string& file1, const std::string& file2) {
    std::ifstream f1(file1, std::ios::binary);
    std::ifstream f2(file2, std::ios::binary);

    if (!f1.is_open() || !f2.is_open()) {
        throw std::runtime_error("Could not open one or both files.");
    }

    // Compare file sizes
    f1.seekg(0, std::ios::end);
    f2.seekg(0, std::ios::end);
    if (f1.tellg() != f2.tellg()) {
        return false; // Different sizes
    }

    f1.seekg(0, std::ios::beg);
    f2.seekg(0, std::ios::beg);

    // Compare byte by byte
    char byte1, byte2;
    while (f1.get(byte1) && f2.get(byte2)) {
        if (byte1 != byte2) {
            return false; // Files differ
        }
    }

    return true; // Files are identical
}

TEST_CASE("Test Monte Carlo Run", "[test_simulation][MonteCarloRun]") {
    // Reference configuration
    std::string config_path = std::string(PROJECT_ROOT_DIR) + "/tests/config/initconf.xyz";
    configuration cfg = ReadTrimCFG(config_path);

    // Params 
    double T;
    int tau;
    int cycles;
    int tw;
    double p_flip;
    std::vector<std::string> observables = {"U", "MSD", "Fs"};
    std::string out;
    int n_log;
    int n_lin;
    std::string params_path = std::string(PROJECT_ROOT_DIR) + "/tests/params/params.json";

    ReadJSONParams(
        params_path, out, N, T, tau, tw, cycles, n_log, n_lin, p_flip
    );

    int seed = 12345; // Specific seed for reproducibility
    srand(seed);

    // Making the out directory
    MakeOutDir(out, params_path);

    // Running
    MonteCarloRun(cfg, T, tau, cycles, tw, p_flip, observables, out, n_log, n_lin);

    SECTION("Check if the dynamics is correct") {
        std::string last_cfg = out + "configs/cfg_500.xy";
        std::string ref_cfg = std::string(PROJECT_ROOT_DIR) + "/tests/reference/cfg_100.xy";
        REQUIRE(AreFilesIdentical(last_cfg, ref_cfg)==true);
    }
    
    SECTION("Check if the observables evolution is correct") {
        std::string obs = out + "obs.txt";
        std::string ref_obs = std::string(PROJECT_ROOT_DIR) + "/tests/reference/obs.txt";
        REQUIRE(AreFilesIdentical(obs, ref_obs)==true);
    }

    // Cleanup: Remove the output directory
    fs::remove_all(out);
}