#include "swap.h"
std::vector <configuration> cfgsCycles;

// Monte Carlo Simulation loop
void MC(std::string out, int n_log, int n_lin){
    int dataCounter=0;
    int cycle;
    int cycleCounter = 0;
    
    // Building snapshots list (log-spaced)
    std::vector < std::pair <double, double>> pairs;
    std::vector <double> samplePoints, twPoints;
    double endingPoints[cycles], linPoints[n_lin];
    double exponents = log10(tau)/(n_log-1);

    for(int c=0; c<cycles; c++){
        for (int x = 0; x < n_log; x++){
            double value = tw*c + floor(pow(10,exponents*(x)));
            std::pair <double,double> p = {value, c};
            int f = std::count(pairs.begin(), pairs.end(), p);
            if(f==0){
                pairs.emplace_back(value, c);
            // this if condition is actually relevent because of the floor function
            }
        }
    }

    // Sorting
    std::sort(pairs.begin(), pairs.end());
    for (auto p: pairs){
        samplePoints.push_back(p.first); twPoints.push_back(p.second);
    }

    // Ending points
    for(int c=0;c<cycles;c++){
        endingPoints[c] = c*tw + tau;
    }
    // Linspaced points
    for (int k=1;k<=n_lin;k++){
        linPoints[k] = (tau/(n_lin))*k;
    }

    // File writing
    std::ofstream log_obs, log_cfg; 
    std::string out_cfg = out + "configs/";
    log_obs.open(out + "obs.txt");
    log_obs << "t" << " " << "cycle";
    for (std::string obs: allObs){
        log_obs << " " << obs;
    } log_obs << std::endl;
    log_obs << std::scientific << std::setprecision(8);
    // creating configs dir
    fs::create_directory(out_cfg); 

    for(int t = 1; t <= steps; t++){
        // Checking whether to update the neighbours list
        cfg.CheckNL();
    
        // Updating reference observables
        if((t-1)%tw == 0 && cycleCounter < cycles){
            cfgsCycles.push_back(cfg); cycleCounter++;
        } 

        // // Writing observables to text file
        int lin = std::count(linPoints, linPoints+n_lin, 1.0*t);
        int log = std::count(samplePoints.begin(), samplePoints.end(), 1.0*t);

        if(lin>0){ // checking if linear saving time
            // Configs
            log_cfg.open(out_cfg + "cfg_" + std::to_string(t) + ".xy");
            log_cfg << std::scientific << std::setprecision(8);
            for (int i = 0; i<N; i++){
                log_cfg << cfg.S[i] << " " << cfg.Xfull[i] << " " << cfg.Yfull[i] << " " << cfg.Zfull[i] << std::endl;
            }
            log_cfg.close();
        }

        if(log>0){ // checking if log saving time
            for(int s=0; s<log; s++){
                // looping different eventual tws
                cycle = twPoints[dataCounter];
                // Configs
                if(! fs::exists (out_cfg + "cfg_" + std::to_string(t) + ".xy")){
                    log_cfg.open(out_cfg + "cfg_" + std::to_string(t) + ".xy");
                    log_cfg << std::scientific << std::setprecision(8);
                    for (int i = 0; i<N; i++){
                        log_cfg << cfg.S[i] << " " << cfg.Xfull[i] << " " << cfg.Yfull[i] << " " << cfg.Zfull[i] << std::endl;
                    }
                    log_cfg.close();
                } 
                // observables
                log_obs << t << " " << cycle;
                for (std::string obs: allObs){
                    log_obs << " " << whichObs(obs, cycle);
                } log_obs << std::endl;

                dataCounter++;
            }  
        };
        // Doing the MC
        for (int i = 0; i < N; i++){
            if (ranf() > p_flip) TryDisp(floor(ranf()*N)); //Displacement probability 0.8
            else TryFlip(floor(ranf()*N)); //Flip probability 0.2
        }
        
        if((t-1)%100==0) std::cout << (t-1) << std::endl;; // Counting steps
    };
    log_obs.close();
}

//  Tries displacing one particle j by vector dr = (dx, dy, dz)
void TryDisp(int j){
    double dx = (ranf()-0.5)*deltaMax;
    double dy = (ranf()-0.5)*deltaMax;
    double dz = (ranf()-0.5)*deltaMax;
    double Xold = cfg.X[j], Xnew = Pshift(cfg.X[j]+dx); 
    double Yold = cfg.Y[j], Ynew = Pshift(cfg.Y[j]+dy);
    double Zold = cfg.Z[j], Znew = Pshift(cfg.Z[j]+dz);
    // Energy before the displacement
    double V_old = V(j);
    // Energy after the displacement
    cfg.X[j] = Xnew; cfg.Y[j] = Ynew; cfg.Z[j] = Znew;
    cfg.Xfull[j] += dx; cfg.Yfull[j] += dy; cfg.Zfull[j] += dz;
    double V_new = V(j);

    double deltaE = V_new - V_old;
    if (deltaE < 0){
        // pass 
    }
    else if (exp(-deltaE/T) < ranf()){
        cfg.X[j] = Xold; cfg.Y[j] = Yold; cfg.Z[j] = Zold;
        cfg.Xfull[j] -= dx; cfg.Yfull[j] -= dy; cfg.Zfull[j] -= dz;
    }
}

//  Tries swapping two particles diameters in the molecule containing particle j
void TryFlip(int j){
    int a = rand() % 2; int k = BN[j][a]; 
    // Energy of the two clusters before the move attempt
    double V_old = V(X[j],Y[j],Z[j],S[j],j)+V(X[k],Y[k],Z[k],S[k],k);
    // Temporarily saving old configurations
    double Sj_old = S[j]; double Sk_old = S[k];
    S[j] = Sk_old; S[k] = Sj_old;
    // Energy of the two clusters after the move attempt
    double V_new = V(X[j],Y[j],Z[j],S[j],j)+V(X[k],Y[k],Z[k],S[k],k);

    double deltaE = V_new - V_old;
    if (deltaE < 0){
        // pass
    }
    else if (exp(-deltaE/T) < ranf()){
        S[j] = Sj_old; S[k] = Sk_old;
    }
}

// Converts string to real observable
double whichObs(std::string obs, int cycl){
    if (obs=="MSD") return MSD(cfgsCycles[cycl]);
    else if (obs=="U") return VTotal()/(2*N);
    else if (obs=="Fs") return FS(cfgsCycles[cycl]);
}