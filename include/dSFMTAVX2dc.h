#pragma once
#ifndef DSFMTAVX2DC_H
#define DSFMTAVX2DC_H

#include "devavxprng.h"

class options {
public:
    int mexp;
    bool verbose;
    bool fixed;
    int fixedSL1;
    uint64_t seed;
    std::string filename;
    long count;
};

bool parse_opt(options& opt, int argc, char **argv);

int search(options& opt, int count);

#endif
