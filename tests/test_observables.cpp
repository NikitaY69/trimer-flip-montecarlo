#include <catch2/catch.hpp>
#include <cmath>
#include <vector>
#include "globals.hpp"
#include "particles.hpp"
#include "observables.hpp"
#include "utils.hpp"

// Reference configuration
std::string config_path = std::string(PROJECT_ROOT_DIR) + "/tests/config/initconf.xyz";
configuration cfg = ReadTrimCFG(config_path);

TEST_CASE("Test WCAPair function", "[test_observables][WCAPair]") {
    
    SECTION("Check that dumb-particles very close to the cutoff give either 0 or non-zero energy") {
        double epsilon = 0.01; // Small perturbation to the cut-off radius
        double s1 = GENERATE(from_range(std::begin(diameters), std::end(diameters)));
        double s2 = GENERATE(from_range(std::begin(diameters), std::end(diameters)));
        double rc = pow(2., 1./6.)*(s1+s2)/2;

        double result_noninteracting = WCAPair(0.0, 0.0, 0.0, s1, 0.0, 0.0, rc + epsilon, s2);
        REQUIRE(result_noninteracting == 0.0);

        double result_interacting = WCAPair(0.0, rc - epsilon, 0.0, s1, 0.0, 0.0, 0.0, s2);
        REQUIRE(result_interacting != 0.0);
    }

    SECTION("Check that real particles (loaded) give the expected energy") {
        int i1 = 135, i2 = 137, i3 = 453, i4 = 456;
        double s1, s2, s3, s4;
        s1 = diameters[int(cfg.S[i1]-1)]; s2 = diameters[int(cfg.S[i2]-1)];
        s3 = diameters[int(cfg.S[i3]-1)]; s4 = diameters[int(cfg.S[i4]-1)];

        REQUIRE(
            WCAPair(
                cfg.X[i1], cfg.Y[i1], cfg.Z[i1], s1, 
                cfg.X[i2], cfg.Y[i2], cfg.Z[i2], s2
            ) == Approx(6.0807844631057675));

        REQUIRE(
            WCAPair(
                cfg.X[i3], cfg.Y[i3], cfg.Z[i3], s3, 
                cfg.X[i4], cfg.Y[i4], cfg.Z[i4], s4
            ) == Approx(0.0));
    
    }
}

// TEST_CASE("Test FENEPair function", "[FENEPair]") {
//     double result = FENEPair(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0);
//     REQUIRE(result == Approx(expected_value)); // Replace expected_value with the correct value
// }

// TEST_CASE("Test V function", "[V]") {
//     std::string input_path = "path/to/test/configuration.cfg";
//     configuration cfg = load_test_configuration(input_path);
//     double result = V(cfg, 0);
//     REQUIRE(result == Approx(expected_value)); // Replace expected_value with the correct value
// }