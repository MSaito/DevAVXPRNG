/**
 * @file SFMTAVX2dc.cpp
 *
 */

#include "devavxprng.h"
#include "SFMTAVX2search.hpp"
#include "SFMTAVXdc.hpp"

int main(int argc, char** argv) {
    using namespace MTToolBox;
    DCOptions opt(607);
    opt.useSR1 = true;
    opt.fixedSL1 = 19;
    opt.fixedSR1 = 7;
    opt.fixedPerm = 1;
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        return -1;
    }
    return sfmtavx_search<w256_t, SFMTAVX2, 256>(opt, opt.count);
}
