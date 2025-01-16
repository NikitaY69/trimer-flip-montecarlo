cmake_minimum_required(VERSION 3.30)

project(TrimerFlipMC 
        VERSION 1.0
        LANGUAGES CXX)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directory
include_directories("include")

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

FetchContent_MakeAvailable(json indicators)

# Source files
file(GLOB SRC_FILES "src/*.cpp")

# Add executable
add_executable(TrimerFlipMC ${SRC_FILES})

# Link libraries
target_link_libraries(TrimerFlipMC stdc++fs boost_program_options 
                      nlohmann_json::nlohmann_json indicators::indicators)

# Add compiler options
target_compile_options(TrimerFlipMC PRIVATE -O3 -mfma -mbmi2 -flto)