#include "devavxprng.h"
#include "SFMTAVX2search.hpp"
#include "EQOptions.hpp"
#include "SFMTAVXeq.hpp"

using namespace MTToolBox;

int main(int argc, char * argv[])
{
    EQOptions<SFMTAVX2_param> opt;
    if (!opt.parse(argc, argv)) {
        return -1;
    }
    int r
        = sfmtavx_equidistribution<w256_t, SFMTAVX2, SFMTAVX2_param, 256>(opt);
    return r;
}
