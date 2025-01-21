#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <cmath>
#include "globals.hpp"

// Globals that are normally defined in the main file
int N = 3000;
double Size = pow(N/density, 1/3.);