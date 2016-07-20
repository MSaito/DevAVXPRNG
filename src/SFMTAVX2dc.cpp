/**
 * @file SFMTAVX2dc.cpp
 *
 */

#include "devavxprng.h"
#include "SFMTAVX2dc.hpp"

int main(int argc, char** argv) {
    using namespace MTToolBox;
    DCOptions opt(607);
    bool parse = opt.parse(argc, argv);
    if (!parse) {
        return -1;
    }
    return MTToolBox::search(opt, opt.count);
}
