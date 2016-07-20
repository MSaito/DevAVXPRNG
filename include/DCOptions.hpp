#pragma once
#ifndef DCOPTIONS_HPP
#define DCOPTIONS_HPP
/**
 * @file DCOptions.hpp
 */
#include "devavxprng.h"
#include <string>
#include <sstream>
#include <fstream>

namespace MTToolBox {
    /**
     * Dynamic creator options
     */
    class DCOptions {
    public:
        int mexp;
        bool verbose;
        bool fixed;
        uint64_t seed;
        long count;
        int min_mexp;

        DCOptions(int min_mexp) {
            mexp = 0;
            verbose = false;
            fixed = false;
            seed = (uint64_t)clock();
            count = 1;
            this->min_mexp = min_mexp;
        }
#if 1
        void d_p() {
            using namespace std;
            cout << "mexp:" << dec << mexp << endl;
            cout << "verbose:" << verbose << endl;
            cout << "fixed:" << fixed << endl;
            cout << "seed:" << dec << seed << endl;
            cout << "count:" << dec << count << endl;
            cout << "min_mexp:" << dec << min_mexp << endl;
        }
#endif
        /**
         * command line option parser
         * @param opt a structure to keep the result of parsing
         * @param argc number of command line arguments
         * @param argv command line arguments
         * @param start default start value
         * @return command line options have error, or not
         */
        bool parse(int argc, char **argv) {
            using namespace std;
            int c;
            bool error = false;
            string pgm = argv[0];
            static struct option longopts[] = {
                {"verbose", no_argument, NULL, 'v'},
                {"fixed", no_argument, NULL, 'x'},
                {"count", required_argument, NULL, 'c'},
                {"seed", required_argument, NULL, 's'},
                {NULL, 0, NULL, 0}};
            errno = 0;
            for (;;) {
                c = getopt_long(argc, argv, "vxs:c:", longopts, NULL);
                if (error) {
                    break;
                }
                if (c == -1) {
                    break;
                }
                switch (c) {
                case 's':
                    seed = strtoull(optarg, NULL, 0);
                    if (errno) {
                        error = true;
                        cerr << "seed must be a number" << endl;
                    }
                    break;
                case 'v':
                    verbose = true;
                    break;
                case 'x':
                    fixed = true;
                    break;
                case 'c':
                    count = strtoll(optarg, NULL, 10);
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
                const int allowed_mexp[] = {607, 1279, 2281, 3217,
                                            4253, 4423,
                                            9689, 9941, 11213, 19937,
                                            21701, 23209, 44497, -1};
                if (! errno) {
                    if (mexp < min_mexp) {
                        error = true;
                    }
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
                        if (allowed_mexp[i] >= min_mexp) {
                            cerr << dec << allowed_mexp[i] << " ";
                        }
                    }
                    cerr << endl;
                }
                this->mexp = mexp;
            }
            if (error) {
                output_help(pgm);
                return false;
            }
            return true;
        }
    private:
        /**
         * showing help message
         * @param pgm program name
         */
        void output_help(std::string& pgm) {
            using namespace std;
            cerr << "usage:" << endl;
            cerr << pgm
                 << " [-s seed] [-v] [-c count]"
                 << " [-x]"
                 << " mexp"
                 << endl;
            static string help_string1 = "\n"
                "--verbose, -v        Verbose mode. Output parameters,"
                " calculation time, etc.\n"
                "--count, -c count    Output count. The number of parameters"
                " to be outputted.\n"
                "--seed, -s seed      seed of randomness.\n"
                "--fixed, -x          use prefixed parameter.\n"
                "mexp                 mersenne exponent.\n"
                ;
            cerr << help_string1 << endl;
        }
    };
}
#endif // DCOPTIONS_HPP
