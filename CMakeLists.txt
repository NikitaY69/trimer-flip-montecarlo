cmake_minimum_required(VERSION 3.30)

project(TFMC 
        VERSION 1.0
        LANGUAGES CXX)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define the project root directory
add_definitions(-DPROJECT_ROOT_DIR=\"${CMAKE_SOURCE_DIR}\")
set(PROJECT_ROOT_DIR "${CMAKE_SOURCE_DIR}")

# Include directories
include_directories("include")

include(FetchContent)

# Fetch the indicators library
FetchContent_Declare(
    indicators
    GIT_REPOSITORY https://github.com/p-ranav/indicators.git
    GIT_TAG v2.3
)

# Fetch Catch2
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.7
)

FetchContent_MakeAvailable(indicators Catch2)

# Find Boost
set(Boost_USE_STATIC_LIBS ON) # Use static libraries
set(Boost_USE_MULTITHREADED ON) # Use multithreaded libraries
set(Boost_USE_STATIC_RUNTIME OFF) # Use dynamic runtime

find_package(Boost REQUIRED COMPONENTS program_options json filesystem)

# Find Boost
# set(Boost_NO_BOOST_CMAKE ON)
# find_package(Boost REQUIRED COMPONENTS program_options)
# if(Boost_FOUND)
#     include_directories(${Boost_INCLUDE_DIRS})
#     link_directories(${Boost_LIBRARY_DIRS})
#     message(STATUS "Boost found: ${Boost_INCLUDE_DIRS}")
# else()
#     message(FATAL_ERROR "Boost not found")
# endif()

# Ensure the GCC LTO-friendly archiver and ranlib are used
set(CMAKE_AR "gcc-ar")
set(CMAKE_RANLIB "gcc-ranlib")

# Source files for the library
file(GLOB SRC_FILES "src/particles.cpp" "src/observables.cpp" "src/simulation.cpp"
                    "src/utils.cpp")

# Create a static library
add_library(TFMC_lib STATIC ${SRC_FILES})

# Link libraries to the static library
target_link_libraries(TFMC_lib Boost::program_options Boost::json Boost::filesystem indicators::indicators)

# Add compiler options for the library
target_compile_options(TFMC_lib PUBLIC -O3 -mfma -mbmi2 -flto)

# Add executable for the main application
add_executable(TFMC src/main.cpp)
target_link_libraries(TFMC TFMC_lib)

# Install the executable
install(TARGETS TFMC 
        RUNTIME DESTINATION bin)

# Enable testing
enable_testing()

# Configure the params file
configure_file(
    ${CMAKE_SOURCE_DIR}/tests/params/params_template.json
    ${CMAKE_SOURCE_DIR}/tests/params/params.json
    @ONLY
)

# Add tests
file(GLOB TEST_FILES "tests/*.cpp")

add_executable(TFMC_tests ${TEST_FILES})
target_link_libraries(TFMC_tests TFMC_lib Catch2::Catch2)

# Register tests with ctest
add_test(NAME test_utils COMMAND TFMC_tests [test_utils] -r compact)
add_test(NAME test_particles COMMAND TFMC_tests [test_particles] -r compact)
add_test(NAME test_observables COMMAND TFMC_tests [test_observables] -r compact)
add_test(NAME test_simulation COMMAND TFMC_tests [test_simulation] -r compact)