#pragma once
#ifndef DSFMTAVXEQMPI_HPP
#define DSFMTAVXEQMPI_HPP
/**
 * @file dSFMTAVXeqmpi.hpp
 */

#include "devavxprng.h"
#include "AlgorithmDSFMTEquidistribution.hpp"
//#include <MTToolBox/AlgorithmEquidistribution.hpp>
//#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
//#include <MTToolBox/period.hpp>
//#include <NTL/GF2X.h>
#include "Annihilate.hpp"
#include "EQOptions.hpp"

namespace MTToolBox {

    template<typename U, typename G, typename P, int bitWidth>
    int dsfmtavxmpi_equidist(EQOptions<P>& opt,
                             int rank, int num_process)
    {
        using namespace std;
        G sf(opt.params);
        U wseed;
        wseed.u64[0] = opt.seed;
        sf.seed(wseed);
        Annihilate<G, U, uint64_t> annihilate;
        if (!annihilate.anni(sf)) {
            return -1;
        }

        DSFMTInfo info;
        info.bitSize = bitWidth; // IMPORTANT
        info.elementNo = bitWidth / 64;
        for (int v = 1; v <= 64; v++) {
            if (v % num_process == rank) {
                int veq = calc_dSFMT_equidist<U, G>(v, sf, info,
                                                    opt.params.mexp);
                int d = opt.params.mexp / v - veq;
                cout << "k(" << dec << v << ") = " << dec << veq;
                cout << "\td(" << dec << v << ") = " << dec << d << endl;
            }
        }
        return 0;
    }
}

#endif // DSFMTAVXEQMPI_HPP
