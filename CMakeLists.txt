cmake_minimum_required(VERSION 3.7.2)

project(TrimerFlipMC 
        VERSION 1.0
        LANGUAGES CXX)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directory
include_directories("include")

# Source files
file(GLOB SRC_FILES "src/*.cpp")

# Add executable
add_executable(TrimerFlipMC ${SRC_FILES})

# Link libraries
target_link_libraries(TrimerFlipMC stdc++fs boost_program_options)

# Add compiler options
target_compile_options(TrimerFlipMC PRIVATE -O3 -mfma -mbmi2 -flto)