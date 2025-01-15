#include <vector>
#include <fstream>
#include <sstream>
#include "globals.hpp"
#include "io.hpp"

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

// Read polydisperse configs
configuration ReadPolyCFG(std::string input){
    configuration C;
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
            C.S[i] = cfg[i][0]; 
            C.X[i] = Pshift(cfg[i][1]); C.Y[i] = Pshift(cfg[i][2]); C.Z[i] = Pshift(cfg[i][3]);
            C.X0[i] = C.X[i]; C.Xfull[i] = C.X[i];
            C.Y0[i] = C.Y[i]; C.Yfull[i] = C.Y[i];
            C.Z0[i] = C.Z[i]; C.Zfull[i] = C.Z[i];
            i++;}
        input_file.close();
        
    } else {
        std::string error = input + " not found";
        throw error;
    }
}