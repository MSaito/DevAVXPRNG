/**
 * @file SFMTAVX2dc.cpp
 *
 */

#include "devavxprng.h"
#include <string>
#include <sstream>
#include <fstream>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <NTL/GF2X.h>
#include "SFMTAVX2search.hpp"
#include "AlgorithmSIMDEquidistribution.hpp"
#include "Annihilate.hpp"
#include "SFMTAVX2dc.h"

using namespace std;
using namespace MTToolBox;
using namespace NTL;

/**
 * search parameters using all_in_one function in the file search_all.hpp
 * @param opt command line options
 * @param count number of parameters user requested
 * @return 0 if this ends normally
 */
int search(options& opt, int count) {
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

static void output_help(string& pgm);

/**
 * command line option parser
 * @param opt a structure to keep the result of parsing
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param start default start value
 * @return command line options have error, or not
 */
bool parse_opt(options& opt, int argc, char **argv) {
    opt.verbose = false;
    opt.mexp = 0;
    opt.count = 1;
    opt.fixed = false;
    opt.seed = (uint64_t)clock();
    opt.filename = "";
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"file", required_argument, NULL, 'f'},
        {"count", required_argument, NULL, 'c'},
        {"fixed", no_argument, NULL, 'x'},
        {"seed", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vxs:f:c:", longopts, NULL);
        if (error) {
            break;
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 's':
            opt.seed = strtoull(optarg, NULL, 0);
            if (errno) {
                error = true;
                cerr << "seed must be a number" << endl;
            }
            break;
        case 'v':
            opt.verbose = true;
            break;
        case 'x':
            opt.fixed = true;
            break;
        case 'f':
            opt.filename = optarg;
            break;
        case 'c':
            opt.count = strtoll(optarg, NULL, 10);
            if (errno) {
                error = true;
                cerr << "count must be a number" << endl;
            }
            break;
        case '?':
        default:
            error = true;
            break;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc < 1) {
        error = true;
    } else {
        long mexp = strtol(argv[0], NULL, 10);
        static const int allowed_mexp[] = {607, 1279, 2281, 3217, 4253, 4423,
                                           9689, 9941, 11213, 19937,
                                           21701, 23209, 44497, -1};
        if (! errno) {
            bool found = false;
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                if (mexp == allowed_mexp[i]) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                error = true;
            }
        }
        if (errno || error){
            error = true;
            cerr << "mexp must be one of ";
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                cerr << dec << allowed_mexp[i] << " ";
            }
            cerr << endl;
        }
        opt.mexp = mexp;
    }
    if (!opt.filename.empty()) {
        ofstream ofs(opt.filename.c_str());
        if (ofs) {
            ofs.close();
        } else {
            error = true;
            cerr << "can't open file:" << opt.filename << endl;
        }
    }
    if (error) {
        output_help(pgm);
        return false;
    }
    return true;
}

/**
 * showing help message
 * @param pgm program name
 */
static void output_help(string& pgm) {
    cerr << "usage:" << endl;
    cerr << pgm
         << " [-s seed] [-v] [-c count]"
         << " [-f outputfile]"
         << " mexp"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output parameters, calculation time, etc.\n"
"--file, -f filename  Parameters are outputted to this file. without this\n"
"                     option, parameters are outputted to standard output.\n"
"--count, -c count    Output count. The number of parameters to be outputted.\n"
"--seed, -s seed      seed of randomness.\n"
"mexp                 mersenne exponent.\n"
        ;
    cerr << help_string1 << endl;
}

#if !defined(NO_MAIN)
int main(int argc, char** argv) {
    options opt;
    bool parse = parse_opt(opt, argc, argv);
    if (!parse) {
        return -1;
    }
    return search(opt, opt.count);
}
#endif
