#include "swap.h"
double x_max; // Maximum value of r/s for real neighbours

// Method to calculate center of mass coordinates
double configuration::GetCM_coord(int coord_index){
    double a = 0;
    std::vector <double>* coordinates;
    (coord_index == 0) ? coordinates = &Xfull : 
    (coord_index == 1) ? coordinates = &Yfull : 
    (coord_index == 2) ? coordinates = &Zfull : 
    throw std::invalid_argument("Invalid coord_index");
    for (int i=0; i<N; i++){
        a += (*coordinates)[i];
    } return a/N;
}

// Method to calculate the verlet lists associated to each particle
void configuration::UpdateNL(){
    NL.clear(); NL = std::vector < std::vector <int> > (N);
    for (int j=0; j<N-1; j++){
        for (int i=j+1; i<N; i++){
            double xij = bcs(X[i], X[j]); 
            double yij = bcs(Y[i], Y[j]); 
            double zij = bcs(Z[i], Z[j]);
            double rij2 = (xij*xij)+(yij*yij)+(zij*zij);
            if (rij2 < rNL && i != j){
                NL[j].push_back(i);
                NL[i].push_back(j);
            }
        }
    }
}

// Retrieves bonded particles for all particles (done only once)
// ATM it is implicit that configurations are written in the trimers index order
void configuration::GetBonds(){
    BN.clear(); BN = std::vector < std::vector<int> > (N);
    for (int i=0; i<N; i+=3){
        BN[i].push_back(i+1); BN[i].push_back(i+2);
        BN[i+1].push_back(i); BN[i+1].push_back(i+2);
        BN[i+2].push_back(i); BN[i+2].push_back(i+1);
    }
}
    
//  Calculates difference of a and b while applying periodic boundary conditions
double bcs(double a, double b) {return Size/2 - std::abs(std::abs(a-b)-Size/2);}

// Shifts coordinate inside main box
double Pshift(double a){
    return fmod(a, Size);//- Size*floor((a+Size/2)/Size);
}

// // Computes the nearest neighbours list
// void UpdateNN(int t0){
//     NN.clear(); NN = std::vector < std::vector <int> > (N);
//     if (t0 == 0){
//         x_max = 1.485;
//     }
//     else{
//         x_max = 1.7;
//     }
//     for (int j=0; j<N-1; j++){
//         for (int i=j+1; i<N; i++){
//             double sigmaij = (S[i]+S[j])*(1-0.2*std::abs(S[i]-S[j]))/2;
//             double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]); double zij = bcs(Z[i], Z[j]);
//             double rij = sqrt((xij*xij)+(yij*yij)+(zij*zij));
//             if (rij < x_max*sigmaij && i != j){
//                 NN[j].push_back(i);
//                 NN[i].push_back(j);
//             }
//         } 
//     }
// }