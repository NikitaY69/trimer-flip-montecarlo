#ifndef IO_H
#define IO_H

#include <string>
#include "particles.hpp"

configuration ReadTrimCFG(std::string input);

void WriteTrimCFG(const configuration& cfg, std::string output);

#endif