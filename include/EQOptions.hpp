#pragma once
#ifndef EQOPTIONS_HPP
#define EQOPTIONS_HPP
/**
 * @file EQOptions.hpp
 */
#include "devavxprng.h"
#include <stdlib.h>
#include <getopt.h>

namespace MTToolBox {
    template<typename P>
    class EQOptions {
    public:
        bool verbose;
        uint64_t seed;
        P params;

        EQOptions() {
            using namespace std;
            verbose = false;
            seed = (uint64_t)clock();
        }

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
            errno = 0;
            string pgm = argv[0];
            static struct option longopts[] = {
                {"verbose", no_argument, NULL, 'v'},
                {"seed", required_argument, NULL, 's'},
                {NULL, 0, NULL, 0}};
            for (;;) {
                c = getopt_long(argc, argv, "vs:", longopts, NULL);
                if (error) {
                    break;
                }
                if (c == -1) {
                    break;
                }
                switch (c) {
                case 'v':
                    verbose = true;
                    break;
                case 's':
                    seed = strtoull(optarg, NULL, 0);
                    if (errno) {
                        error = true;
                        cerr << "seed must be a number" << endl;
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
                params.readFromString(argv[0]);
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
                     << " [-v] [-s] \""
                     << params.get_header()
                     << "\""
                     << endl;
                cerr << "\n"
                     << "--verbose, -v        Verbose mode. Output detailed "
                     << "information.\n"
                     << "--seed, -s seed      seed of randomness.\n";
        }
    };
}
#endif // EQOPTIONS_HPP
