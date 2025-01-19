cmake_minimum_required(VERSION 3.10)

project(TFMC 
        VERSION 1.0
        LANGUAGES CXX)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories("include" "tests")

include(FetchContent)

# Fetch json library
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.2
)

# Fetch the indicators library
FetchContent_Declare(
    indicators
    GIT_REPOSITORY https://github.com/p-ranav/indicators.git
    GIT_TAG master  # You can specify a specific commit or branch if necessary
)

# Fetch Catch2
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.7 # Specify the version you want to use
)

FetchContent_MakeAvailable(json indicators Catch2)

# Source files
file(GLOB SRC_FILES "src/particles.cpp" "src/observables.cpp" "src/simulation.cpp"
                    "src/utils.cpp")

# Create a static library for the common source files
add_library(TFMC_lib STATIC ${SRC_FILES})

# Link libraries to the static library
target_link_libraries(TFMC_lib stdc++fs boost_program_options 
                      nlohmann_json::nlohmann_json indicators::indicators)

# Add executable for the main application
add_executable(TFMC src/main.cpp)
target_link_libraries(TFMC TFMC_lib)

# Add compiler options
target_compile_options(TFMC PRIVATE -O3 -mfma -mbmi2 -flto)

# Install the executable
install(TARGETS TFMC 
        RUNTIME DESTINATION bin)

# Enable testing
enable_testing()

# Add tests
file(GLOB TEST_FILES "tests/test_observables.cpp")
add_executable(TFMC_tests ${TEST_FILES})
target_link_libraries(TFMC_tests TFMC_lib Catch2::Catch2)
add_test(NAME TFMC_tests COMMAND TFMC_tests)

# Define the project root directory
add_definitions(-DPROJECT_ROOT_DIR=\"${CMAKE_SOURCE_DIR}\")