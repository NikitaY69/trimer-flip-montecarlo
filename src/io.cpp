#include "swap.h"

// Read polydisperse configs
void ReadPolyCFG(std::string input){
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
            S[i] = cfg[i][0]; 
            X[i] = Pshift(cfg[i][1]); Y[i] = Pshift(cfg[i][2]); Z[i] = Pshift(cfg[i][3]);
            X0[i] = X[i]; Xfull[i] = X[i]; Xref[i] = X[i]; 
            Y0[i] = Y[i]; Yfull[i] = Y[i]; Yref[i] = Y[i];
            Z0[i] = Z[i]; Zfull[i] = Z[i]; Zref[i] = Z[i];
            Sref[i] = S[i];
            i++;}
        input_file.close();
        
    } else {
        std::string error = input + " not found";
        throw error;
    }
}

// Read trimer configs
void ReadTrimCFG(std::string input){
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
            mol_index[i] = cfg[i][0];
            type = cfg[i][1]-1; 
            S[i] = diameters[type];
            X[i] = Pshift(cfg[i][2]); Y[i] = Pshift(cfg[i][3]); Z[i] = Pshift(cfg[i][4]);
            X0[i] = X[i]; Xfull[i] = X[i]; Xref[i] = X[i]; 
            Y0[i] = Y[i]; Yfull[i] = Y[i]; Yref[i] = Y[i];
            Z0[i] = Z[i]; Zfull[i] = Z[i]; Zref[i] = Z[i];
            Sref[i] = S[i];
            i++;}
        input_file.close();
        
    } else {
        std::string error = input + " not found";
        throw error;
    }
}