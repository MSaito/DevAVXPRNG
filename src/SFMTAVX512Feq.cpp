#include "devavxprng.h"
#include "SFMTAVX512Fsearch.hpp"
#include "EQOptions.hpp"
#include "SFMTAVXeq.hpp"

using namespace MTToolBox;

int main(int argc, char * argv[])
{
    EQOptions<SFMTAVX512F_param> opt;
    if (!opt.parse(argc, argv)) {
        return -1;
    }
    int r
        = sfmtavx_equidistribution<w512_t, SFMTAVX512F,
                                   SFMTAVX512F_param, 512>(opt);
    return r;
}
