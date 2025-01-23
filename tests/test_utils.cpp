#include <catch2/catch.hpp>
#include <string>
#include <vector>
#include "globals.hpp"
#include "utils.hpp"

TEST_CASE("Test ParseCMDLine function", "[test_utils][ParseCMDLine]") {
    const char* argv[] = {
        "program", 
        "--init", "input.cfg", 
        "--params", "params.json", 
        "--observables", "MSD", "FS"};
    int argc = sizeof(argv) / sizeof(argv[0]);

    std::string input;
    std::string params;
    std::vector<std::string> observables;
    int seed;

    bool result = ParseCMDLine(argc, argv, input, params, observables);
    REQUIRE(result == true);
    REQUIRE(input == "input.cfg");
    REQUIRE(params == "params.json");
    REQUIRE(observables.size() == 2);
}