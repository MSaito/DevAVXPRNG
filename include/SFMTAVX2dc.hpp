#pragma once
#ifndef SFMTAVX2DC_HPP
#define SFMTAVX2DC_HPP
/**
 * @file SFMTAVX2dc.cpp
 *
 */

#include "devavxprng.h"
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <NTL/GF2X.h>
#include "SFMTAVX2search.hpp"
#include "AlgorithmSIMDEquidistribution.hpp"
#include "Annihilate.hpp"
#include "DCOptions.hpp"

namespace MTToolBox {
    /**
     * search parameters using all_in_one function in the file search_all.hpp
     * @param opt command line options
     * @param count number of parameters user requested
     * @return 0 if this ends normally
     */
    int search(DCOptions& opt, int count) {
        using namespace std;
        using namespace NTL;

        MersenneTwister mt(opt.seed);
        SFMTAVX2 g(opt.mexp);
        if (opt.fixed) {
            g.setFixed(true);
        }

        cout << "seed = " << dec << opt.seed << endl;
        if (opt.verbose) {
            time_t t = time(NULL);
            cout << "search start at " << ctime(&t);
        }
        AlgorithmReducibleRecursionSearch<w256_t, uint32_t> ars(g, mt);
        int i = 0;
        AlgorithmCalculateParity<w256_t, SFMTAVX2> cp;
        Annihilate<SFMTAVX2, w256_t, uint32_t> annihilate;
        cout << "# " << g.getHeaderString() << ", delta32, delta64"
             << endl;
        while (i < count) {
            if (ars.start(opt.mexp * 100)) {
                GF2X irreducible = ars.getIrreducibleFactor();
                GF2X characteristic = ars.getCharacteristicPolynomial();
                GF2X quotient = characteristic / irreducible;
                if (deg(irreducible) != opt.mexp) {
                    cout << "error" << endl;
                    return -1;
                }
                cp.searchParity(g, irreducible);
#if 1
                w256_t seed = {{1, 0, 0, 0, 0, 0, 0, 0}};
                g.seed(seed);
                if (!annihilate.anni(g)) {
                    return -1;
                }
                int delta32 = 0;
                int delta64 = 0;
                int veq32[32];
                SIMDInfo info;
                info.bitMode = 32;
                info.elementNo = 8;
                info.bitSize = 256;
                info.fastMode = true;
                delta32
                    = calc_SIMD_equidistribution<w256_t, SFMTAVX2>
                    (g, veq32, 32, info, opt.mexp);
                int veq64[64];
                info.bitMode = 64;
                info.elementNo = 4;
                delta64
                    = calc_SIMD_equidistribution<w256_t, SFMTAVX2>
                    (g, veq64, 64, info, opt.mexp);
#endif
                cout << g.getParamString();
                cout << dec << delta32 << "," << delta64;
                cout << "," << dec << veq64[63] << endl;
                i++;
            } else {
                cout << "search failed" << endl;
                break;
            }
        }
        if (opt.verbose) {
            time_t t = time(NULL);
            cout << "search end at " << ctime(&t) << endl;
        }
        return 0;
    }

}
#endif // SFMTAVX2DC_HPP
