#include "swap.h"

double dXCM, dYCM, dZCM;
int dataCounter=0;
int cycle;

// Monte Carlo Simulation loop
void MC(std::string out, int n_log, int n_lin){
    int cycleCounter = 0;
    double deltaX[N], deltaY[N], deltaZ[N], deltaR2[N], R2Max = 0;
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
        // Updating NL
        if((t-1) % 1 == 0) {//Change number?
            // every 150 steps we check if we need to update the NL
            for (int i = 0; i < N; i++){
                deltaX[i] = bcs(X[i],X0[i]);
                deltaY[i] = bcs(Y[i],Y0[i]);
                deltaZ[i] = bcs(Z[i],Z0[i]);
                deltaR2[i] = deltaX[i]*deltaX[i] + deltaY[i]*deltaY[i] + deltaZ[i]*deltaZ[i];
            R2Max = std::max_element(deltaR2,deltaR2+N)[0];
            }
            if(R2Max > RUpdate){
                // std::cout << (t-1) << std::endl;
                UpdateNL();
                R2Max = 0;
                for(int j = 0; j < N; j++){
                    X0[j] = X[j];
                    Y0[j] = Y[j];
                    Z0[j] = Z[j];
                }
            }
        }
    
        // Updating reference observables
        if((t-1)%tw == 0 && cycleCounter < cycles){
            UpdateAge(cycleCounter); cycleCounter++;
        } 

        // // Writing observables to text file
        int lin = std::count(linPoints, linPoints+n_lin, 1.0*t);
        int log = std::count(samplePoints.begin(), samplePoints.end(), 1.0*t);

        if(lin>0){ // checking if linear saving time
            // Configs
            log_cfg.open(out_cfg + "cfg_" + std::to_string(t) + ".xy");
            log_cfg << std::scientific << std::setprecision(8);
            for (int i = 0; i<N; i++){
                log_cfg << S[i] << " " << Xfull[i] << " " << Yfull[i] << " " << Zfull[i] << std::endl;
            }
            log_cfg.close();
        }

        if(log>0){ // checking if log saving time
            // UpdateNN(t); // updating nearest neighbours
            dXCM = 0; dYCM = 0, dZCM = 0;
            for (int i=0;i<N;i++){
                double dX = Xfull[i]-Xref[i], dY = Yfull[i]-Yref[i], dZ = Zfull[i]-Zref[i];
                dXCM += dX; dYCM += dY; dZCM += dZ;
            } dXCM /= N; dYCM /= N; dZCM /= N;

            for(int s=0; s<log; s++){
                // looping different eventual tws
                cycle = twPoints[dataCounter];
                // Configs
                if(! fs::exists (out_cfg + "cfg_" + std::to_string(t) + ".xy")){
                    log_cfg.open(out_cfg + "cfg_" + std::to_string(t) + ".xy");
                    log_cfg << std::scientific << std::setprecision(8);
                    for (int i = 0; i<N; i++){
                        log_cfg << S[i] << " " << Xfull[i] << " " << Yfull[i] << " " << Zfull[i] << std::endl;
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
        
        // if((t-1)%100==0) std::cout << (t-1) << std::endl;; // Counting steps
    };
    log_obs.close();
}

//  Tries displacing one particle j by vector dr = (dx, dy, dz)
void TryDisp(int j){
    double dx = (ranf()-0.5)*deltaMax;
    double dy = (ranf()-0.5)*deltaMax;
    double dz = (ranf()-0.5)*deltaMax;
    double Xnew = Pshift(X[j]+dx);
    double Ynew = Pshift(Y[j]+dy);
    double Znew = Pshift(Z[j]+dz);
    double deltaE = V(Xnew, Ynew, Znew, S[j], j) - V(X[j], Y[j], Z[j], S[j], j);
    // why is the modulus function not in deltaE ?
    if (deltaE < 0){
        // Xnew = fmod(X[j],Size);
        X[j] = Xnew; //Check modulus function
        Y[j] = Ynew;
        Z[j] = Znew;
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
        Zfull[j] = Zfull[j]+dz;
    }
    else if (exp(-deltaE/T) > ranf()){
        X[j] = Xnew;
        Y[j] = Ynew;
        Z[j] = Znew;
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
        Zfull[j] = Zfull[j]+dz;
    }
}

//  Tries swapping two particles diameters in the molecule containing particle j
void TryFlip(int j){
    int a = rand() % 2; int k = BN[j][a]; 
    double deltaE = V(X[j],Y[j],Z[j],S[k],j)+V(X[k],Y[k],Z[k],S[j],k)-V(X[j],Y[j],Z[j],S[j],j)-V(X[k],Y[k],Z[k],S[k],k);
    if (deltaE < 0){
        double Rnew = S[k];
        S[k] = S[j];
        S[j] = Rnew;
        // swapCount[j] += 1; swapCount[k] += 1;
    }
    else if (exp(-deltaE/T) > ranf()){
        double Rnew = S[k];
        S[k] = S[j];
        S[j] = Rnew;
        // swapCount[j] += 1; swapCount[k] += 1;
    }
}

// Converts string to real observable
double whichObs(std::string obs, int cycl){
    if (obs=="MSD") return MSD();
    else if (obs=="U") return VTotal()/(2*N);
    else if (obs=="Cb") return CB(cycl);
    else if (obs=="Fs") return FS(cycl);
}