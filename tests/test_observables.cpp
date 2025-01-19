#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <cmath>
#include <vector>
#include "globals.hpp"
#include "observables.hpp"
#include "test_helpers.hpp"

// Reference configuration
// cfg = LoadTestConfiguration();

TEST_CASE("Test WCAPair function", "[WCAPair]") {
    double epsilon = 0.01; // Small perturbation to the cut-off radius
    double s1 = GENERATE(from_range(std::begin(diameters), std::end(diameters)));
    double s2 = GENERATE(from_range(std::begin(diameters), std::end(diameters)));
    
    SECTION("Check that dumb-particles very close to the cutoff give either 0 or non-zero energy") {
        double rc = pow(2., 1./6.)*(s1+s2)/2;

        double result_noninteracting = WCAPair(0.0, 0.0, 0.0, s1, 0.0, 0.0, rc + epsilon, s2);
        REQUIRE(result_noninteracting == 0.0);

        double result_interacting = WCAPair(0.0, rc - epsilon, 0.0, s1, 0.0, 0.0, 0.0, s2);
        REQUIRE(result_interacting != 0.0);
    }
}

//     SECTION("Check that real particles (loaded) give the expected energy") {
//         int i1 = 0, i2 = 1, s1, s2;
//         sj = diameters[int(cfg.S[i1]-1)]; sk = diameters[int(cfg.S[i2]-1)]
//         double result = WCAPair(
//             cfg.X[particle1], cfg.Y[particle1], cfg.Z[particle1], cfg.S[particle1], 
//             cfg.X[particle2], cfg.Y[particle2], cfg.Z[particle2], cfg.S[particle2]
//         );
//         REQUIRE(result == Approx(expected_value)); // Replace expected_value with the correct value
//     }
    
// }

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