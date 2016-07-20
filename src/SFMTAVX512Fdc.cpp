/**
 * @file SFMTAVX512Fdc.cpp
 */
#include "devavxprng.h"
#include "DCOptions.hpp"
#include "SFMTAVX512Fdc.hpp"

using namespace MTToolBox;

int main(int argc, char** argv) {
    DCOptions opt(1279);
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        return -1;
    }
    return search(opt, opt.count);
}
