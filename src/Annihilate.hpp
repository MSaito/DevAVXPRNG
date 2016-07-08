#pragma once
#ifndef ANNIHILATE_HPP
#define ANNIHILATE_HPP

#include "devavxprng.h"
#include <NTL/GF2X.h>
#include <MTToolBox/period.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include "w256.hpp"

namespace MTToolBox {
    /**
     * @tparam G generator
     * @tparam U simd vector like w256_t
     * @tparam V base rng width like uint32_t
     */
    template<typename G, typename U, typename V>
    class Annihilate {
    public:
        bool anni(G& sf) {
            using namespace NTL;
            using namespace std;
            //G gen(sf);
            GF2X poly;
            minpoly(poly, sf);
            GF2X irreducible = poly;
            if (!hasFactorOfDegree(irreducible, sf.getMexp())) {
                cout << "error does not have factor of degree(0) "
                     << dec << sf.getMexp()
                     << endl;
                cout << "deg poly = " << dec << deg(poly) << endl;
                cout << "deg irreducible = " << dec << deg(irreducible) << endl;
                cout << sf.getParamString() << endl;
                sf.d_p();
                cout << "poly:" << endl;
                cout << poly << endl;
                cout << "irreducible:" << endl;
                cout << irreducible << endl;
                return false;
            }
#if 1
            calcCharacteristicPolynomial(&sf, poly);
            if (deg(poly) != sf.bitSize()) {
                getLCMPoly(poly, sf);
            }
#endif
            int degp = deg(poly);
            //printBinary(stdout, poly);
            GF2X quotient = poly / irreducible;
#if defined(DEBUG)
            cout << "deg irreducible = " << dec << deg(irreducible) << endl;
            cout << "deg characteristic = " << dec << deg(poly)
                 << endl;
            cout << "deg quotient = " << dec << deg(quotient) << endl;
#endif
            GF2X tmp = quotient * irreducible;
            if (tmp != poly) {
                cout << "quotient * irreducible != poly" << endl;
                return false;
            }
            annihilate<U, V>(&sf, quotient);
            minpoly(poly, sf);
            if (poly != irreducible) {
                cout << "annihilate failed" << endl;
                cout << "deg poly = " << dec << degp << endl;
                cout << "deg irreducible = " << dec << deg(irreducible) << endl;
                cout << "deg quotient = " << dec << deg(quotient) << endl;
                cout << "after annihilate deg minpoly = " << dec << deg(poly)
                     << endl;
                return false;
            }
#if defined(DEBUG)
            cout << "after annihilate deg poly = " << dec << deg(poly) << endl;
#endif
            if (!hasFactorOfDegree(poly, sf.getMexp())) {
                cout << "error does not have factor of degree(1) "
                     << dec << sf.getMexp() << endl;
                cout << "deg poly = " << dec << deg(poly) << endl;
                return false;
            }
            return true;
        }
        void getLCMPoly(NTL::GF2X& lcm, const G& sf) {
            using namespace NTL;
            using namespace std;
            G gen(sf);
            lungLCM(lcm, gen);
            int bitSize = gen.bitSize();
            GF2X poly;
            for (int i = 0; i < bitSize; i++) {
                gen.setOneBit(i);
#if 0
                for (int j = 0; j < 256; j++) {
                    minPolyLung(poly, gen, j);
                    //minpoly(poly, gen, j);
                    LCM(lcm, lcm, poly);
                    if (deg(lcm) == bitSize) {
                        return;
                    }
                }
#else
                minPolyLung(poly, gen, 0);
                //minpoly(poly, gen);
                LCM(lcm, lcm, poly);
                if (deg(lcm) == bitSize) {
                    return;
                }
#endif
            }
        }
    private:
        void lungLCM(NTL::GF2X& lcm, G& sf) {
            using namespace NTL;
            using namespace std;
            GF2X poly;
            for (int i = 0; i < 256; i++) {
                minPolyLung(poly, sf, i);
                LCM(lcm, lcm, poly);
            }
        }

        void minPolyLung(NTL::GF2X& poly, G& sf, int pos) {
            using namespace NTL;
            using namespace std;
            Vec<GF2> v;
            int size = sf.bitSize();
            v.SetLength(2 * size);
            for (int i = 0; i < 2 * size; i++) {
                sf.generate();
                w256_t w = sf.getParityValue();
                v[i] = getBitOfPos(w, pos);
            }
            MinPolySeq(poly, v, size);
        }
    };
}
#endif // ANNIHILATE_HPP
