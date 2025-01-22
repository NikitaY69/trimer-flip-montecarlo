#include <cmath> // For math operations
#include "globals.hpp"
#include "observables.hpp"

// Universal constants
const double pi = 3.14159265358979323846;

//  Calculates the pairwise WCA potential between two particles
double WCAPair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2){
    // int idx1 = s1-1, idx2 = s2-1;
    // double sigmaij = (diameters[idx1]+diameters[idx2])/2;
    double sigmaij = (s1+s2)/2;
    double sigma2 = sigmaij*sigmaij;
    double rc2 = pow(2., 1./3.) * sigma2;
    double xij = MinimumImageDistance(x1, x2); 
    double yij = MinimumImageDistance(y1, y2); 
    double zij = MinimumImageDistance(z1, z2);
    double rij2 = (xij*xij) + (yij*yij) + (zij*zij);
    if (rij2 > rc2) return 0;
    else {
        double a2 = sigma2/rij2; double a6 = a2*a2*a2;
        return 4*(a6*a6-a6+0.25);
    }
}

//  Calculates the pairwise FENE potential between two particles
double FENEPair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2){
    // int idx1 = s1-1, idx2 = s2-1;
    // double sigmaij = (diameters[idx1]+diameters[idx2])/2;
    double sigmaij = (s1+s2)/2;
    double sigma2 = sigmaij*sigmaij;
    double kij = 30/sigma2;
    double R02 = 1.5*1.5*sigma2;
    double xij = MinimumImageDistance(x1, x2); 
    double yij = MinimumImageDistance(y1, y2); 
    double zij = MinimumImageDistance(z1, z2);
    double rij2 = (xij*xij) + (yij*yij) + (zij*zij);
    if (rij2 > R02) return 0;
    else {
        return -0.5*kij*R02*log(1-rij2/R02);
    }
}

//  Calculates potential associated to particle j
double V(const configuration& cfg, int j){
    double total = 0, sj, sk;
    for (int k: cfg.neighbours_list[j]){
        sj = diameters[int(cfg.S[j]-1)]; sk = diameters[int(cfg.S[k]-1)];
        total += WCAPair(cfg.X[j], cfg.Y[j], cfg.Z[j], sj, 
                         cfg.X[k], cfg.Y[k], cfg.Z[k], sk);
    }
    for (int k: cfg.bonded_neighbours[j]){
        sj = diameters[int(cfg.S[j]-1)]; sk = diameters[int(cfg.S[k]-1)];
        total += FENEPair(cfg.X[j], cfg.Y[j], cfg.Z[j], sj, 
                          cfg.X[k], cfg.Y[k], cfg.Z[k], sk);
    } return total;
}

//  Calculates total system energy (double)
double VTotal(const configuration& cfg){
    double vTot = 0;
    for (int j = 0; j < N; j++)
        vTot += V(cfg, j);
    return vTot;
}

//  Calculates avg. mean square displacements
double MSD(const configuration& cfg, const configuration& cfg0){
    double sum = 0, deltaX, deltaY, deltaZ;
        for (int i = 0; i < N; i++){
            deltaX = cfg.Xfull[i]-cfg0.Xfull[i]; deltaX -= (cfg.XCM-cfg0.XCM);
            deltaY = cfg.Yfull[i]-cfg0.Yfull[i]; deltaY -= (cfg.YCM-cfg0.YCM);
            deltaZ = cfg.Zfull[i]-cfg0.Zfull[i]; deltaZ -= (cfg.ZCM-cfg0.ZCM);
            sum += deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
    }
    return sum/N;
}

// Correlation functions

//  Calculates the intermediate self-scattering function
double FS(const configuration& cfg, const configuration& cfg0){
    double dotProduct;
    double q = 2*pi/sigmaMax;
    double sum = 0, deltaX, deltaY, deltaZ;
    int ang = 360;
    
    for (int theta=0; theta<ang; theta++){
        for (int i = 0; i < N; i++){
            deltaX = cfg.Xfull[i]-cfg0.Xfull[i]; deltaX -= (cfg.XCM-cfg0.XCM);
            deltaY = cfg.Yfull[i]-cfg0.Yfull[i]; deltaY -= (cfg.YCM-cfg0.YCM);
            deltaZ = cfg.Zfull[i]-cfg0.Zfull[i]; deltaZ -= (cfg.ZCM-cfg0.ZCM);
            dotProduct = q*((cos(theta*pi/180)*deltaX)+(sin(theta*pi/180)*deltaY));
            sum += cos(dotProduct);
        }
    }
    return sum/(ang*N);
}