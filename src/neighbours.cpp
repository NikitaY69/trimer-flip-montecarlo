#include "swap.h"
double r_step = 35.0/(nr-1);
double x_max; // Maximum value of r/s for real neighbours

//  Calculates difference of a and b while applying periodic boundary conditions
double bcs(double a, double b) {return Size/2 - std::abs(std::abs(a-b)-Size/2);}

// Shifts coordinate inside main box
double Pshift(double a){
    return a - Size*floor((a+Size/2)/Size);
}

// Retrieves bonded particles for all particles (done only once)
// ATM it is implicit that configurations are written in the trimers index order
std::vector < std::vector<int> > GetBonds(){
    std::vector < std::vector<int> > vec(N);
    for (int i=0; i<N; i+=3){
        vec[i].push_back(i+1); vec[i].push_back(i+2);
        vec[i+1].push_back(i); vec[i+1].push_back(i+2);
        vec[i+2].push_back(i); vec[i+2].push_back(i+1);
    }
    return vec;
}

// Computes the pseudo-interacting neighbours list
void UpdateNL(){
    NL.clear(); NL = std::vector < std::vector <int> > (N);
    for (int j=0; j<N-1; j++){
        for (int i=j+1; i<N; i++){
            double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]); double zij = bcs(Z[i], Z[j]);
            double rij2 = (xij*xij)+(yij*yij)+(zij*zij);
            if (rij2 < rNL && i != j){
                NL[j].push_back(i);
                NL[i].push_back(j);
            }
        }
    }
}

// Computes the nearest neighbours list
void UpdateNN(int t0){
    NN.clear(); NN = std::vector < std::vector <int> > (N);
    if (t0 == 0){
        x_max = 1.485;
    }
    else{
        x_max = 1.7;
    }
    for (int j=0; j<N-1; j++){
        for (int i=j+1; i<N; i++){
            double sigmaij = (S[i]+S[j])*(1-0.2*std::abs(S[i]-S[j]))/2;
            double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]); double zij = bcs(Z[i], Z[j]);
            double rij = sqrt((xij*xij)+(yij*yij)+(zij*zij));
            if (rij < x_max*sigmaij && i != j){
                NN[j].push_back(i);
                NN[i].push_back(j);
            }
        } 
    }
}

// Computes the per-radius neighbours list
void UpdateRL(){
    RL.clear(); RL = std::vector < std::vector < std::vector <int>>> 
    (N, std::vector < std::vector <int>>(nr));
    for (int i=0;i<N-1;i++){
        for (int j=i+1; j<N; j++){
            double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]); double zij = bcs(Z[i], Z[j]);
            double rij = sqrt((xij*xij)+(yij*yij)+(zij*zij));
            for (int k=0; k<nr; k++){
                double r = k*r_step;
                if (rij<=r){
                    RL[j][k].push_back(i);
                    RL[i][k].push_back(j);
                }
            }
        }
    }
}
