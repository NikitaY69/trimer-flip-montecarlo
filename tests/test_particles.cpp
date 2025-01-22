#include <catch2/catch.hpp>
#include "globals.hpp"
#include "particles.hpp"
#include "utils.hpp"

// Reference configuration
// std::string config_path = std::string(PROJECT_ROOT_DIR) + "/tests/config/initconf.xyz";
// configuration cfg = ReadTrimCFG(config_path);

// Test the bcs function
TEST_CASE("Test MinimumImageDistance function", "[test_particles][MinimumImageDistance]") {
    double epsilon = 0.01;
    REQUIRE(MinimumImageDistance(0., Size/2+epsilon) == Approx(Size/2-epsilon));
    REQUIRE(MinimumImageDistance(Size/2-epsilon, Size/2+epsilon) == Approx(2*epsilon));
    REQUIRE(MinimumImageDistance(0., Size-epsilon) == Approx(epsilon));
}

// Test the Pshift function
TEST_CASE("Test ShiftInMainBox function", "[test_particles][ShiftInMainBox]") {
    double epsilon = 0.01;

    SECTION("Check Pshift at boundaries") {
        REQUIRE(ShiftInMainBox(0.0) == Approx(0.0));
        REQUIRE(ShiftInMainBox(Size) == Approx(0.0));
        REQUIRE(ShiftInMainBox(-Size) == Approx(0.0));
    }

    SECTION("Check that Pshift is periodic") {
        REQUIRE(ShiftInMainBox(Size+epsilon) == Approx(epsilon));
        REQUIRE(ShiftInMainBox(-epsilon) == Approx(Size-epsilon));
    }
    
}

// Test the UpdateCM_coord method
// TEST_CASE("Test UpdateCM_coord method", "[test_particles][UpdateCM_coord]") {
//     configuration cfg;
//     cfg.N = 3;
//     cfg.Xfull = {1.0, 2.0, 3.0};
//     cfg.Yfull = {4.0, 5.0, 6.0};
//     cfg.Zfull = {7.0, 8.0, 9.0};

//     cfg.UpdateCM_coord();

//     REQUIRE(cfg.XCM == Approx(2.0));
//     REQUIRE(cfg.YCM == Approx(5.0));
//     REQUIRE(cfg.ZCM == Approx(8.0));
// }

// // Test the UpdateNL method
// TEST_CASE("Test UpdateNL method", "[test_particles][UpdateNL]") {
//     configuration cfg;
//     cfg.N = 3;
//     cfg.X = {0.0, 1.0, 2.0};
//     cfg.Y = {0.0, 1.0, 2.0};
//     cfg.Z = {0.0, 1.0, 2.0};
//     cfg.X0 = cfg.X;
//     cfg.Y0 = cfg.Y;
//     cfg.Z0 = cfg.Z;

//     cfg.UpdateNL();

//     REQUIRE(cfg.NL.size() == 3);
//     REQUIRE(cfg.NL[0].size() == 2);
//     REQUIRE(cfg.NL[1].size() == 2);
//     REQUIRE(cfg.NL[2].size() == 2);
// }

// // Test the CheckNL method
// TEST_CASE("Test CheckNL method", "[CheckNL]") {
//     configuration cfg;
//     cfg.N = 3;
//     cfg.X = {0.0, 1.0, 2.0};
//     cfg.Y = {0.0, 1.0, 2.0};
//     cfg.Z = {0.0, 1.0, 2.0};
//     cfg.X0 = {0.0, 0.5, 1.5};
//     cfg.Y0 = {0.0, 0.5, 1.5};
//     cfg.Z0 = {0.0, 0.5, 1.5};

//     cfg.CheckNL();

//     REQUIRE(cfg.NL.size() == 3);
//     REQUIRE(cfg.NL[0].size() == 2);
//     REQUIRE(cfg.NL[1].size() == 2);
//     REQUIRE(cfg.NL[2].size() == 2);
// }