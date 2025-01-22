#include <cmath>
#include <algorithm>
#include "globals.hpp"
#include "particles.hpp"

// Neighbours lists parameters
const double r_cutoff = pow(2., 1./6.) * sigmaMax; // Cutoff radius for calculating potential
const double r_skin = 0.7; // Margin added to the neighbours list radius
const double neighbours_radius_squared = pow(r_cutoff+r_skin,2); // NL radius squared
const double maximum_displacement_before_update_squared = pow(r_skin,2)/4; // When R2Max exceeds this, update NL

// Method to calculate center of mass coordinates
void configuration::UpdateCM_coord(){
    XCM = 0; YCM = 0; ZCM = 0;
    // std::vector <double>* coordinates;
    // (coord_index == 0) ? coordinates = &Xfull : 
    // (coord_index == 1) ? coordinates = &Yfull : 
    // (coord_index == 2) ? coordinates = &Zfull : 
    // throw std::invalid_argument("Invalid coord_index");
    for (int i=0; i<N; i++){
        XCM += Xfull[i]; YCM += Yfull[i]; ZCM += Zfull[i];
    } XCM /= N; YCM /= N; ZCM /= N;
}

// Method to calculate the verlet lists associated to each particle
void configuration::UpdateNL(){
    neighbours_list.clear(); neighbours_list = std::vector < std::vector <int> > (N);
    for (int j=0; j<N-1; j++){
        for (int i=j+1; i<N; i++){
            double xij = bcs(X[i], X[j]); 
            double yij = bcs(Y[i], Y[j]); 
            double zij = bcs(Z[i], Z[j]);
            double rij2 = (xij*xij)+(yij*yij)+(zij*zij);
            if (rij2 < neighbours_radius_squared && i != j){
                neighbours_list[j].push_back(i);
                neighbours_list[i].push_back(j);
            }
        }
    }
}

// Retrieves bonded particles for all particles (done only once)
// ATM it is implicit that configurations are written in the trimers index order
void configuration::GetBonds(){
    bonded_neighbours.clear(); bonded_neighbours = std::vector < std::vector<int> > (N);
    for (int i=0; i<N; i+=3){
        bonded_neighbours[i].push_back(i+1); bonded_neighbours[i].push_back(i+2);
        bonded_neighbours[i+1].push_back(i); bonded_neighbours[i+1].push_back(i+2);
        bonded_neighbours[i+2].push_back(i); bonded_neighbours[i+2].push_back(i+1);
    }
}

// Checking whether to update the neighbours list
void configuration::CheckNL(){
    double deltaX, deltaY, deltaZ, maximum_displacement = 0;
    std::vector <double> deltaR2(N);
    for (int i = 0; i < N; i++){
        deltaX = bcs(X[i],X0[i]);
        deltaY = bcs(Y[i],Y0[i]);
        deltaZ = bcs(Z[i],Z0[i]);
        deltaR2[i] = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
    } maximum_displacement = *(std::max_element(deltaR2.begin(), deltaR2.end()));
    if(maximum_displacement > maximum_displacement_before_update_squared){
        // std::cout << (t-1) << std::endl;
        UpdateNL();
        maximum_displacement = 0;
        X0 = X; Y0 = Y; Z0 = Z;
    }
}
//  Calculates difference of a and b while applying periodic boundary conditions
double bcs(double a, double b) {return Size/2 - std::abs(std::abs(a-b)-Size/2);}

// Shifts coordinate inside main box
double Pshift(double a){
    double shift = fmod(a, Size);//- Size*floor((a+Size/2)/Size);
    if (shift < 0) {shift += Size;}
    return shift;
}