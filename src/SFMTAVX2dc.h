#pragma once
#ifndef SFMTAVX2DC_H
#define SFMTAVX2DC_H

#include "devavxprng.h"

class options {
public:
    int mexp;
    bool verbose;
    bool fixed;
    uint64_t seed;
    std::string filename;
    long count;
};

bool parse_opt(options& opt, int argc, char **argv);

int search(options& opt, int count);

#endif // SFMTAVX2DC_H
