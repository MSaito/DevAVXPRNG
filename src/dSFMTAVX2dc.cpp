/**
 * @file dSFMTAVX2dc.cpp
 *
 * @brief The main function of parameter generator of dSFMTAVX2.
 *
 */
#include "devavxprng.h"
#include "dSFMTAVX2search.hpp"
#include "DCOptions.hpp"
#include "dSFMTAVXdc.hpp"

int main(int argc, char** argv) {
    using namespace MTToolBox;
    DCOptions opt(607);
    opt.fixed = true;
    opt.fixedValue = 19;
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        return -1;
    }
    return search<w256_t, dSFMTAVX2, 256>(opt, opt.count);
}
