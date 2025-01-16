#include <cmath>
#include <fstream>
#include <experimental/filesystem>
#include <algorithm>
#include <iostream>
#include <indicators/progress_bar.hpp>
#include "globals.hpp"
#include "simulation.hpp"
#include "utils.hpp"
#include "observables.hpp"

namespace fs = std::experimental::filesystem;

// Constants
const double deltaMax = 0.12; // Max particle displacement

// Progress bar
indicators::ProgressBar bar{
    indicators::option::BarWidth{50},
    indicators::option::Start{"["},
    indicators::option::Fill{"="},
    indicators::option::Lead{">"},
    indicators::option::Remainder{" "},
    indicators::option::End{"]"},
    indicators::option::ForegroundColor{indicators::Color::cyan},
    indicators::option::ShowPercentage{true},
    indicators::option::ShowElapsedTime{true},
};

// Monte Carlo Simulation loop
void MonteCarloRun(configuration& cfg, double T, int tau, int cycles, int tw, double p_flip, 
        std::vector <std::string>& observables, std::string& out, int n_log, int n_lin){
            
    double steps = tw*(cycles-1)+tau;
    int dataCounter=0;
    int cycle;
    int cycleCounter = 0;
    std::vector <configuration> cfgsCycles;
    configuration* cfg0;

    // Building snapshots list
    std::vector <int> logpoints, twpoints, linpoints;

    // Logspaced
    std::vector < std::pair <int, int>> log_and_tws = GetLogspacedSnapshots(cycles, tau, tw, n_log);
    for (auto p: log_and_tws){
        logpoints.push_back(p.first); twpoints.push_back(p.second);
    }
    // Linspaced
    linpoints = GetLinspacedSnapshots(cycles, tau, tw, n_lin);

    // File writing
    std::string out_cfg = out + "configs/";
    std::ofstream log_obs = MakeObsFile(observables, out + "obs.txt");
    // creating configs dir
    fs::create_directory(out_cfg); 

    // First neighbours
    cfg.GetBonds(); cfg.UpdateNL();

    // Monte Carlo sweeps
    for(int t = 1; t <= steps; t++){
        // Checking whether to update the neighbours list
        cfg.CheckNL();
    
        // Updating reference observables
        if((t-1)%tw == 0 && cycleCounter < cycles){
            cfg.UpdateCM_coord();
            cfgsCycles.push_back(cfg); cycleCounter++;
        } 

        // // Writing observables to text file
        int lin = std::count(linpoints.begin(), linpoints.end(), t);
        int log = std::count(logpoints.begin(), logpoints.end(), t);

        if(lin>0){ // checking if linear saving time
            // Configs
            if(! fs::exists (out_cfg + "cfg_" + std::to_string(t) + ".xy")){
                WriteTrimCFG(cfg, out_cfg + "cfg_" + std::to_string(t) + ".xy");
            }
        }

        if(log>0){ // checking if log saving time
            cfg.UpdateCM_coord();
            for(int s=0; s<log; s++){
                // looping different eventual tws
                cycle = twpoints[dataCounter];
                cfg0 = &cfgsCycles[cycle];
                // Configs
                if(! fs::exists (out_cfg + "cfg_" + std::to_string(t) + ".xy")){
                    WriteTrimCFG(cfg, out_cfg + "cfg_" + std::to_string(t) + ".xy");
                } 
                // Observables
                WriteObs(cfg, *cfg0, t, cycle, observables, log_obs);

                dataCounter++;
            }  
        };
        // Doing the MC
        for (int i = 0; i < N; i++){
            if (ranf() > p_flip) TryDisp(cfg, floor(ranf()*N), T); //Displacement probability 0.8
            else TryFlip(cfg, floor(ranf()*N), T); //Flip probability 0.2
        }
        
        if((t-1)%(tau/100)==0) bar.tick();
        
    };
    log_obs.close();
}

//  Tries displacing one particle j by vector dr = (dx, dy, dz)
void TryDisp(configuration& cfg, int j, double T){
    double dx = (ranf()-0.5)*deltaMax;
    double dy = (ranf()-0.5)*deltaMax;
    double dz = (ranf()-0.5)*deltaMax;
    double Xold = cfg.X[j], Xnew = Pshift(cfg.X[j]+dx); 
    double Yold = cfg.Y[j], Ynew = Pshift(cfg.Y[j]+dy);
    double Zold = cfg.Z[j], Znew = Pshift(cfg.Z[j]+dz);
    // Energy before the displacement
    double V_old = V(cfg, j);
    // Energy after the displacement
    cfg.X[j] = Xnew; cfg.Y[j] = Ynew; cfg.Z[j] = Znew;
    cfg.Xfull[j] += dx; cfg.Yfull[j] += dy; cfg.Zfull[j] += dz;
    double V_new = V(cfg, j);

    double deltaE = V_new - V_old;
    if (deltaE < 0){
        // pass 
    }
    else if (exp(-deltaE/T) < ranf()){
        cfg.X[j] = Xold; cfg.Y[j] = Yold; cfg.Z[j] = Zold;
        cfg.Xfull[j] -= dx; cfg.Yfull[j] -= dy; cfg.Zfull[j] -= dz;
    }
}

// Observables-only run
void ComputeObservables(int tau, int cycles, int tw,  
        std::vector <std::string>& observables, std::string& out, int n_log){
    
    double steps = tw*(cycles-1)+tau;
    int dataCounter=0;
    int cycle;
    int cycleCounter = 0;
    std::vector <configuration> cfgsCycles;
    configuration cfg;
    configuration* cfg0;

    // Building snapshots list
    std::vector <int> logpoints, twpoints;

    // Logspaced
    std::vector < std::pair <int, int>> log_and_tws = GetLogspacedSnapshots(cycles, tau, tw, n_log);
    for (auto p: log_and_tws){
        logpoints.push_back(p.first); twpoints.push_back(p.second);
    }

    // File writing
    std::string out_cfg = out + "configs/";
    std::ofstream log_obs = MakeObsFile(observables, out + "obs.txt");

    // Looping over the saved snapshots
    for(int t: logpoints){

        cfg = ReadTrimCFG(out_cfg + "cfg_" + std::to_string(t) + ".xy");
        cfg.UpdateNL();
        cfg.UpdateCM_coord();
        // std::cout << t << std::endl;
        // Updating reference observables
        if((t-1)%tw == 0 && cycleCounter < cycles){
            // std::cout << "updating" << std::endl;
            cfgsCycles.push_back(cfg); cycleCounter++;
        } 
        
        cycle = twpoints[dataCounter];
        cfg0 = &cfgsCycles[cycle];

        // Observables
        WriteObs(cfg, *cfg0, t, cycle, observables, log_obs);

        dataCounter++;
        bar.tick();
   
    };
    log_obs.close();
}

//  Tries swapping two particles diameters in the molecule containing particle j
void TryFlip(configuration& cfg, int j, double T){
    int a = rand() % 2; int k = cfg.BN[j][a]; 
    // Energy of the two clusters before the move attempt
    double V_old = V(cfg, j) + V(cfg, k);
    // Temporarily saving old configurations
    double Sj_old = cfg.S[j]; double Sk_old = cfg.S[k];
    cfg.S[j] = Sk_old; cfg.S[k] = Sj_old;
    // Energy of the two clusters after the move attempt
    double V_new = V(cfg, j) + V(cfg, k);

    double deltaE = V_new - V_old;
    if (deltaE < 0){
        // pass
    }
    else if (exp(-deltaE/T) < ranf()){
        cfg.S[j] = Sj_old; cfg.S[k] = Sk_old;
    }
}

