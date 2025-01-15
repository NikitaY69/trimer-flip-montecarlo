#include <sstream>
#include <iomanip>
#include "globals.hpp"
#include "io.hpp"
#include "observables.hpp"

// Read trimer configs
configuration ReadTrimCFG(std::string input){
    configuration C;
    int type;
    std::string line;
    std::ifstream input_file(input);
    if (input_file.is_open()){
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
            type = cfg[i][1]-1; 
            C.S[i] = (diameters[type]);
            C.X[i] = Pshift(cfg[i][2]); C.Y[i] = Pshift(cfg[i][3]); C.Z[i] = Pshift(cfg[i][4]);
            C.X0[i] = C.X[i]; C.Xfull[i] = C.X[i];
            C.Y0[i] = C.Y[i]; C.Yfull[i] = C.Y[i];
            C.Z0[i] = C.Z[i]; C.Zfull[i] = C.Z[i];
            i++;}
        input_file.close();
        return C;

    } else {
        std::string error = input + " not found";
        throw error;
    }
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