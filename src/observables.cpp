#include "swap.h"

//  Calculates the pairwise purely repulsive potential between two particles
double RepulsivePair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2){
    double sigmaij = (s1+s2)*(1-0.2*std::abs(s1-s2))/2;
    double sigma2 = sigmaij*sigmaij;
    double rc2 = 1.25 * 1.25 * sigma2;
    double xij = bcs(x1, x2); double yij = bcs(y1, y2); double zij = bcs(z1, z2);
    double rij2 = (xij*xij) + (yij*yij) + (zij*zij);

    if (rij2 > rc2) return 0;
    else {
        double a2 = rij2/sigma2; double a4 = a2*a2;
        return (1/(a4*a4*a4))+c0+(c2*a2)+(c4*a4);
    }
}

//  Calculates the pairwise WCA potential between two particles
double WCAPair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2){
    double sigmaij = (s1+s2)/2;
    double sigma2 = sigmaij*sigmaij;
    double rc2 = pow(2., 1./3.) * sigma2;
    double xij = bcs(x1, x2); double yij = bcs(y1, y2); double zij = bcs(z1, z2);
    double rij2 = (xij*xij) + (yij*yij) + (zij*zij);
    if (rij2 > rc2) return 0;
    else {
        double a2 = sigma2/rij2; double a6 = a2*a2*a2;
        return 4*(a6*a6-a6+0.25);
    }
}

//  Calculates the pairwise FENE potential between two particles
double FENEPair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2){
    double sigmaij = (s1+s2)/2;
    double sigma2 = sigmaij*sigmaij;
    double k; double R0;
    if ((s1 == 0.9 && s2 == 1.0) || (s1 == 1.0 && s2 == 0.9)){
        k = 33.241; R0 = 1.425;
    } else if ((s1 == 1.0 && s2 == 1.1) || (s1 == 1.1 && s2 == 1.0)){
        k = 27.210884; R0 = 1.575;
    } else if ((s1 == 0.9 && s2 == 1.1) || (s1 == 1.1 && s2 == 0.9)){
        k = 30.0; R0 = 1.5;
    } else{
        std::cout << s1 << " " << s2 << std::endl;
    }
    double kij = k;
    double R02 = R0*R0;
    // double kij = 30/sigma2;
    // double R02 = 1.5*1.5*sigma2;
    double xij = bcs(x1, x2); double yij = bcs(y1, y2); double zij = bcs(z1, z2);
    double rij2 = (xij*xij) + (yij*yij) + (zij*zij);
    if (rij2 > R02) return 0;
    else {
        return -0.5*kij*R02*log(1-rij2/R02);
    }
}

//  Calculates potential associated to particle j
double V(double xj, double yj, double zj, double rj, int j){
    double total = 0;
    for (int k: NL[j]){
        total += WCAPair(xj, yj, zj, rj, X[k], Y[k], Z[k], S[k]);
    }
    for (int k: BN[j]){
        total += FENEPair(xj, yj, zj, rj, X[k], Y[k], Z[k], S[k]);
    } return total;
}

//  Calculates total system energy (double)
double VTotal(){
    double vTot = 0;
    for (int j = 0; j < N; j++)
        vTot += V(X[j], Y[j], Z[j], S[j], j);
    return vTot;
}

//  Calculates avg. mean square displacements
double MSD(){
    double sum = 0, deltaX, deltaY, deltaZ;
        for (int i = 0; i < N; i++){
            deltaX = Xfull[i]-Xref[i];
            deltaY = Yfull[i]-Yref[i];
            deltaZ = Zfull[i]-Zref[i];
            deltaX -= dXCM; deltaY -= dYCM; deltaZ -= dZCM;
            sum += deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
    }
    return sum/N;
}

// Correlation functions

//  Calculates the intermediate self-scattering function
double FS(int cycle){
    double dotProduct;
    double q = 2*pi/sigmaMax;
    double sum = 0, deltaX, deltaY, deltaZ;
    int ang = 360;
    
    for (int theta=0; theta<ang; theta++){
        for (int i = 0; i < N; i++){
            deltaX = Xfull[i]-Xtw[cycle][i];
            deltaY = Yfull[i]-Ytw[cycle][i];
            deltaZ = Zfull[i]-Ztw[cycle][i];
            deltaX -= dXCM; deltaY -= dYCM; deltaZ -= dZCM;
            dotProduct = q*((cos(theta*pi/180)*deltaX)+(sin(theta*pi/180)*deltaY));
            sum += cos(dotProduct);
        }
    }
    return sum/(ang*N);
}

// Computes the bond-breaking correlation function (local)
double CBLoc(int cycle, int j){
    std::vector<int> intersect;
    std::vector<int> nn0 = NN_tw[cycle][j]; // neighbors at t=0
    std::vector<int> nn = NN[j];
    std::set_intersection(nn0.begin(), nn0.end(), nn.begin(), nn.end(),
                     std::back_inserter(intersect));

    if (nn0.size()==0){
        return 0;
    } else { 
        double frac = intersect.size()/nn0.size();
        return frac;
    } 
}

// Computes the bond-breaking correlation function (averaged)
double CB(int cycle){
    double tot = 0;
    for (int j=0; j<N; j++){
        tot += CBLoc(cycle, j);
    } return tot/N;
}

// Computes the diameter auto-correlation function
double C_sigma(){
    double sigma_m = 1.000218223;
    double C = 0, C0 = 0;
    for (int i=0; i<N; i++){
        double deltaS0 = Sref[i]-sigma_m; double deltaS = S[i]-sigma_m;
        C0 += deltaS0*deltaS0; C += deltaS*deltaS0;
    } return C/C0;
}

// Updates the reference points for the correlation functions
void UpdateAge(int cycle){
    UpdateNN(0); NN_tw.push_back(NN);
    Xtw.push_back(std::vector <double>());
    Ytw.push_back(std::vector <double>());
    Ztw.push_back(std::vector <double>());
    for (int i=0; i<N; i++){
        Xtw[cycle].push_back(Xfull[i]);
        Ytw[cycle].push_back(Yfull[i]);
        Ztw[cycle].push_back(Zfull[i]);
    }
}