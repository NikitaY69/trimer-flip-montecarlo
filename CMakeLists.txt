cmake_minimum_required(VERSION 3.30)

project(TFMC 
        VERSION 1.0
        LANGUAGES CXX)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define the project root directory
add_definitions(-DPROJECT_ROOT_DIR=\"${CMAKE_SOURCE_DIR}\")

# Include directories
include_directories("include")

include(FetchContent)

# Fetch Gulrak/filesystem
FetchContent_Declare(
  ghc_filesystem
  GIT_REPOSITORY https://github.com/gulrak/filesystem.git
)

# Fetch json library
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
)

# Fetch the indicators library
FetchContent_Declare(
    indicators
    GIT_REPOSITORY https://github.com/p-ranav/indicators.git
)

# Fetch Catch2
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG master
)

FetchContent_MakeAvailable(ghc_filesystem json indicators Catch2)

# Find Boost
find_package(Boost REQUIRED COMPONENTS program_options)

if(Boost_FOUND)
    message(STATUS "Boost found: ${Boost_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "Boost not found")
endif()

# Source files
file(GLOB SRC_FILES "src/particles.cpp" "src/observables.cpp" "src/simulation.cpp"
                    "src/utils.cpp")

# Create a static library for the common source files
add_library(TFMC_lib STATIC ${SRC_FILES})

# Link libraries to the static library
target_link_libraries(TFMC_lib ghc_filesystem boost_program_options 
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
file(GLOB TEST_FILES "tests/test_main.cpp" "tests/test_utils.cpp" 
                     "tests/test_particles.cpp" "tests/test_observables.cpp"
                     "tests/test_simulation.cpp")

add_executable(TFMC_tests ${TEST_FILES})
target_link_libraries(TFMC_tests TFMC_lib Catch2::Catch2)
target_compile_options(TFMC_tests PRIVATE -O3 -mfma -mbmi2 -flto)


# Register tests with ctest
add_test(NAME test_utils COMMAND TFMC_tests [test_utils] -r compact)
add_test(NAME test_particles COMMAND TFMC_tests [test_particles] -r compact)
add_test(NAME test_observables COMMAND TFMC_tests [test_observables] -r compact)
add_test(NAME test_simulation COMMAND TFMC_tests [test_simulation] -r compact)