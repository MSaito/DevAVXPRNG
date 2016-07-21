#pragma once
#ifndef SFMTAVXEQMPI_HPP
#define SFMTAVXEQMPI_HPP
/**
 * @file SFMTAVXeqmpi.hpp
 */
#include "devavxprng.h"
#include "AlgorithmSIMDEquidistribution.hpp"
//#include <MTToolBox/AlgorithmEquidistribution.hpp>
//#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
//#include <MTToolBox/period.hpp>
//#include <NTL/GF2X.h>
#include "Annihilate.hpp"
#include "EQOptions.hpp"

namespace MTToolBox {

    template<typename U, typename G, typename P, int bitWidth>
    int sfmtavxmpi_equidist(EQOptions<P>& opt, int rank, int num_process)
    {
        using namespace std;
        G sf(opt.params);
        U wseed;
        wseed.u64[0] = opt.seed;
        sf.seed(wseed);
        Annihilate<G, U, uint32_t> annihilate;
        if (!annihilate.anni(sf)) {
            return -1;
        }

        SIMDInfo info;
        info.bitSize = bitWidth;
        info.fastMode = false;
        sf.reset_reverse_bit();
        info.bitMode = 64;
        info.elementNo = bitWidth / 64;
        if (opt.verbose) {
            cout << "64bit dimension of equidistribution at v-bit accuracy k(v)"
                 << endl;
        }
        for (int v = 1; v <= 64; v++) {
            if (v % num_process == rank) {
                int veq = calc_SIMD_equidist<U, G>(v, sf, info,
                                                   opt.params.mexp);
                int d = opt.params.mexp / v - veq;
                if (opt.verbose) {
                    cout << "k(" << dec << v << ") = " << dec << veq;
                    cout << "\td(" << dec << v << ") = " << dec << d << endl;
                }
            }
        }
        info.bitMode = 32;
        info.elementNo = bitWidth / 32;
        if (opt.verbose) {
            cout << "32bit dimension of equidistribution at v-bit accuracy k(v)"
                 << endl;
        }
        for (int v = 1; v <= 32; v++) {
            if (v % num_process == rank) {
                int veq = calc_SIMD_equidist<U, G>
                    (v, sf, info, opt.params.mexp);
                int d = opt.params.mexp / v - veq;
                if (opt.verbose) {
                    cout << "k(" << dec << v << ") = " << dec << veq;
                    cout << "\td(" << dec << v << ") = " << dec << d << endl;
                }
            }
        }
        return 0;
    }
}
#endif // SFMTAVXEQMPI_HPP
