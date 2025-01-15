#ifndef OBS_H
#define OBS_H

#include "globals.hpp"
#include "particles.hpp"

double WCAPair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2),
       FENEPair(double x1, double y1, double z1, double s1, double x2, double y2, double z2, double s2);

double V(const configuration& cfg, int j), VTotal(const configuration& cfg), 
       MSD(const configuration& cfg, const configuration& cfg0), 
       FS(const configuration& cfg, const configuration& cfg0);

#endif