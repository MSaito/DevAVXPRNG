#include "devavxprng.h"
#include "SFMTAVX2search.hpp"
#include "AlgorithmSIMDEquidistribution.hpp"
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/period.hpp>
#include <NTL/GF2X.h>
#include <errno.h>
#include "Annihilate.hpp"

using namespace MTToolBox;
using namespace std;

class options {
public:
    bool verbose;
    bool reverse;
    uint64_t seed;
    SFMTAVX2_param params;
};

static bool parse_opt(options& opt, int argc, char **argv);
static void output_help(string& pgm);
//static void printBinary(FILE *fp, GF2X& poly);

int main(int argc, char * argv[])
{
    options opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    SFMTAVX2 sf(opt.params);
    //cout << sf.getParamString() << endl;
    w256_t wseed;
    wseed.u64[0] = opt.seed;
    sf.seed(wseed);
    Annihilate<SFMTAVX2, w256_t, uint32_t> annihilate;
    if (!annihilate.anni(sf)) {
        return -1;
    }
    //cout << "annihilate end" << endl;
    int delta32 = 0;
    int delta64 = 0;
    bool lsb = false;
    const char * lsb_str = "";
    if (opt.reverse) {
        sf.set_reverse_bit();
        lsb = true;
        cout << "Equidistribution from LSB" << endl;
        lsb_str = " from LSB";
    }
    //AlgorithmEquidistribution<w256_t, uint32_t> re(sf, 256, opt.params.mexp);
    SIMDInfo info;
    info.fastMode = false;
    info.bitSize = 256;

    info.bitMode = 64;
    info.elementNo = 4;
    if (opt.verbose) {
        cout << "64bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
    }
    for (int v = 1; v <= 64; v++) {
        int veq = calc_SIMD_equidistribution<w256_t, SFMTAVX2, uint32_t>
            (v, sf, info, opt.params.mexp, lsb);
        int d = opt.params.mexp / v - veq;
        if (opt.verbose) {
            cout << "k(" << dec << v << ") = " << dec << veq;
            cout << "\td(" << dec << v << ") = " << dec << d << endl;
        }
        delta64 += d;
    }
    //cout << dec << delta64 << endl;
    info.bitMode = 32;
    info.elementNo = 8;
    sf.reset_reverse_bit();
    if (opt.verbose) {
        cout << "32bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
    }
    for (int v = 1; v <= 32; v++) {
        int veq = calc_SIMD_equidistribution<w256_t, SFMTAVX2, uint32_t>
            (v, sf, info, opt.params.mexp, lsb);
        int d = opt.params.mexp / v - veq;
        if (opt.verbose) {
            cout << "k(" << dec << v << ") = " << dec << veq;
            cout << "\td(" << dec << v << ") = " << dec << d << endl;
        }
        delta32 += d;
    }
    cout << sf.getParamString();
    cout << dec << delta32 << "," << delta64 << endl;
    return 0;
}

/**
 * command line option parser
 * @param opt a structure to keep the result of parsing
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param start default start value
 * @return command line options have error, or not
 */
static bool parse_opt(options& opt, int argc, char **argv) {
    opt.verbose = false;
    opt.seed = (uint64_t)clock();
    opt.reverse = false;
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"seed", required_argument, NULL, 's'},
        {"reverse", no_argument, NULL, 'r'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vrs:", longopts, NULL);
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
        case 'r':
            opt.reverse = true;
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
        opt.params.readFromString(argv[0]);
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
static void output_help(string& pgm)
{
    cerr << "usage:" << endl;
    cerr << pgm
         << " [-v] [-r]"
         << " mexp,pos1,sl1,sr1,msk1,parity1"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output detailed information.\n"
"--reverse, -r        Reverse mode. Calculate equidistribution from LSB.\n"
        ;
    cerr << help_string1 << endl;
}
