#include "devavxprng.h"
#include "SFMTAVX512Fsearch.hpp"
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
    SFMTAVX512F_param params;
};

static bool parse_opt(options& opt, int argc, char **argv);
static void output_help(string& pgm);

int main(int argc, char * argv[])
{
    options opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    SFMTAVX512F sf(opt.params);
    cout << sf.getParamString() << endl;
    w512_t wseed;
    wseed.u64[0] = opt.seed;
    sf.seed(wseed);
    Annihilate<SFMTAVX512F, w512_t, uint32_t> annihilate;
    if (!annihilate.anni(sf)) {
        return -1;
    }
    cout << "annihilate end" << endl;
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
    //AlgorithmEquidistribution<w512_t, uint32_t> re(sf, 512, opt.params.mexp);
    SIMDInfo info;
    info.bitSize = 512;
    int veq32[32];
    info.bitMode = 32;
    info.elementNo = 16;
    info.fastMode = true;
    sf.reset_reverse_bit();
    delta32 = calc_SIMD_equidistribution<w512_t, SFMTAVX512F>
        (sf, veq32, 32, info, opt.params.mexp, lsb);
    cout << "delta32 = " << dec << delta32 << endl;
    int veq64[64];
    info.bitMode = 64;
    info.elementNo = 8;
    delta64 = calc_SIMD_equidistribution<w512_t, SFMTAVX512F>
        (sf, veq64, 64, info, opt.params.mexp, lsb);
    cout << sf.getParamString();
    cout << dec << delta32 << "," << delta64 << endl;
    if (opt.verbose) {
        cout << "32bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
        for (int j = 0; j < 32; j++) {
            cout << "k(" << dec << (j + 1) << ") = " << dec << veq32[j];
            cout << "\td(" << dec << (j + 1) << ") = " << dec
                 << (opt.params.mexp / (j + 1) - veq32[j]) << endl;
        }
        cout << "64bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
        for (int j = 0; j < 64; j++) {
            cout << "k(" << dec << (j + 1) << ") = " << dec << veq64[j];
            cout << "\td(" << dec << (j + 1) << ") = " << dec
                 << (opt.params.mexp / (j + 1) - veq64[j]) << endl;
        }
    }
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
         << " mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,"
         << "parity1,parity2,parity3,parity4"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output detailed information.\n"
"--reverse, -r        Reverse mode. Calculate equidistribution from LSB.\n"
        ;
    cerr << help_string1 << endl;
}
