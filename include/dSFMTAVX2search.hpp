#pragma once
#ifndef DSFMTAVX2SEARCH_HPP
#define DSFMTAVX2SEARCH_HPP

/**
 * @file dSFMTAVX2search.hpp
 *
 * see COPYING
 */

#include "devavxprng.h"
#include <sstream>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/MersenneTwister64.hpp>
#include <MTToolBox/util.hpp>
#include "w256.hpp"

/**
 * @namespace dSFMTAVX2
 * name space for dSFMTAVX2
 */
namespace MTToolBox {
    using namespace NTL;
    using namespace std;

    /**
     * @class dSFMTAVX2_param
     * @brief a class keeping parameters of dSFMTAVX2
     *
     * This class keeps parameters of dSFMTAVX2, and has some
     * method for outputting parameters.
     */
    class dSFMTAVX2_param {
    public:
        int mexp;
        int pos1;
        int sl1;
        int perm;
        w256_t msk1;
        w256_t fix1;
        w256_t parity1;

        dSFMTAVX2_param() {
            mexp = 0;
            pos1 = 0;
            sl1 = 0;
            perm = 1;
            setZero(msk1);
            setZero(fix1);
            setZero(parity1);
        }

        dSFMTAVX2_param(const dSFMTAVX2_param& src) {
            mexp = src.mexp;
            pos1 = src.pos1;
            sl1 = src.sl1;
            perm = src.perm;
            msk1 = src.msk1;
            fix1 = src.fix1;
            parity1 = src.parity1;
        }

        /**
         *
         * @return header line of output.
         */
        const string get_header() const {
            return "mexp, pos1, sl1, perm, msk, fix, parity";
        }

        /**
         *
         * @return string of parameters
         */
        const string get_string() const {
            stringstream ss;
            ss << dec << mexp << ",";
            ss << dec << pos1 << ",";
            ss << dec << sl1 << ",";
            ss << dec << perm << ",";
            ss << msk1 << ",";
            ss << fix1 << ",";
            ss << parity1 << ",";
            string s;
            ss >> s;
            return s;
        }

        /**
         * This method is used for DEBUG.
         * @return string of parameters.
         */
        const string get_debug_string() const {
            stringstream ss;
            ss << "mexp:" << dec << mexp << endl;
            ss << "pos1:" << dec << pos1 << endl;
            ss << "sl1:" << dec << sl1 << endl;
            ss << "perm:" << dec << perm << endl;
            ss << "msk:" << msk1 << endl;
            ss << "fix:" << fix1 << endl;
            ss << "parity:" << parity1 << endl;
            string s;
            ss >> s;
            return s;
        }

        void readFromString(char * str) {
            char * para = str;
            int len = strlen(str);
            mexp = strtoul(para, &para, 10);
            para++;
            pos1 = strtoul(para, &para, 10);
            para++;
            sl1 = strtoul(para, &para, 10);
            para++;
            perm = strtoul(para, &para, 10);
            para++;
            for (int i = 3; i >= 0; i--) {
                msk1.u64[i] = strtoull(para, &para, 16);
                para++;
            }
            if (para >= str + len) {
                return;
            }
            for (int i = 3; i >= 0; i--) {
                fix1.u64[i] = strtoull(para, &para, 16);
                para++;
            }
            if (para >= str + len) {
                return;
            }
            for (int i = 3; i >= 0; i--) {
                parity1.u64[i] = strtoull(para, &para, 16);
                para++;
            }
        }
    };

    /**
     * @class dSFMTAVX2
     * @brief DSFMTAVX2 generator class used for dynamic creation
     *
     * This class is one of the main class of DSFMTAVX2 dynamic creator.
     * This class is designed to be called from programs in MTToolBox.
     */
    class dSFMTAVX2 : public ReducibleGenerator<w256_t> {
    public:
        /**
         * Constructor by mexp.
         * @param mexp Mersenne Exponent
         */
        dSFMTAVX2(int mexp) {
#if defined(DEBUG)
            cout << "dSFMTAVX2 constructor start" << endl;
#endif
            size = (mexp - lung_size) / element_size + 1;
            state = new w256_t[size];
            param.mexp = mexp;
            param.pos1 = 0;
            param.sl1 = 0;
            param.perm = 1;
            MTToolBox::setZero(param.msk1);
            MTToolBox::setZero(param.parity1);
            index = 0;
            start_mode = 0;
            weight_mode = max_weight_mode;
            MTToolBox::setZero(previous);
            MTToolBox::setZero(lung);
            prefix = 0;
            fixedSL1 = 0;
            fixedPerm = 0;
        }

        ~dSFMTAVX2() {
            delete[] state;
        }

        /**
         * The copy constructor.
         * @param src The origin of copy.
         */
        dSFMTAVX2(const dSFMTAVX2& src) : param(src.param) {
#if defined(DEBUG)
            cout << "dSFMTAVX2 constructor start" << endl;
#endif
            size = src.size;
            state = new w256_t[size];
            for (int i = 0; i < size; i++) {
                state[i] = src.state[i];
            }
            lung = src.lung;
            index = src.index;
            start_mode = src.start_mode;
            weight_mode = src.weight_mode;
            previous = src.previous;
            prefix = src.prefix;
            fixedSL1 = src.fixedSL1;
            fixedPerm = src.fixedPerm;
        }

        /**
         * Constructor by parameter.
         * @param src_param
         */
        dSFMTAVX2(const dSFMTAVX2_param& src_param) : param(src_param) {
#if defined(DEBUG)
            cout << "dSFMTAVX2 constructor start" << endl;
#endif
            size = (src_param.mexp - lung_size) / element_size + 1;
            state = new w256_t[size];
            for (int i = 0; i < size; i++) {
                MTToolBox::setZero(state[i]);
            }
            index = 0;
            start_mode = 0;
            weight_mode = max_weight_mode;
            MTToolBox::setZero(previous);
            MTToolBox::setZero(lung);
            prefix = 0;
            fixedSL1 = 0;
            fixedPerm = 0;
        }

        EquidistributionCalculatable<w256_t> * clone() const {
            return new dSFMTAVX2(*this);
        }

        /**
         * This method initialize internal state.
         * This initialization is simple.
         * @param seed seed for initialization
         */
        void seed(w256_t seed) {
            uint64_t * pstate = new uint64_t[(size + 1) * 4];
            for (int i = 0; i < (size + 1) * 4; i++) {
                pstate[i] = 0;
            }
            for (int i = 0; i < 4; i++) {
                pstate[i] = seed.u64[i];
            }
            for (int i = 1; i < (size + 1) * 4; i++) {
                pstate[i] ^= i + UINT64_C(6364136223846793005)
                    * (pstate[i - 1] ^ (pstate[i - 1] >> 62));
            }
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 4; j++) {
                    state[i].u64[j] = pstate[i * 4 + j];
                }
            }
            for (int i = 0; i < 4; i++) {
                lung.u64[i] = pstate[size * 4 + i];
            }
            delete[] pstate;
            index = 0;
            MTToolBox::setZero(previous);
            setup_prefix();
        }

        void do_recursion(w256_t *r, w256_t *a, w256_t *b, w256_t *lung) {
            *lung = permutexvar_epi32(*lung, param.perm);
            *lung ^= *b;
            *lung ^= SL64(*a, param.sl1);
            w256_t t;
            t = SR64(*lung, sr1);
            t ^= *a;
            t ^= *lung & param.msk1;
            *r = t;
        }

        /**
         * Important state transition function.
         */
        void next_state() {
            index = (index + 1) % size;
            do_recursion(&state[index],
                         &state[index],
                         &state[(index + param.pos1) % size],
                         &lung);
        }

        /**
         * Important method, generate new random number
         * @return new pseudo random number
         */
        w256_t generate() {
            next_state();
            w256_t r;
            index = index % size;
            int p = (index + size - 1) % size;
            if (start_mode == 0) {
                r = state[index];
            } else {
                for (int i = 0; i < 4 - start_mode; i++) {
                    r.u64[i] = state[p].u64[i + start_mode];
                }
                int j = 0;
                for (int i = 4 - start_mode; i < 4; i++) {
                    r.u64[i] = state[index].u64[j++];
                }
            }
#if defined(DEBUG) && 0
            if (start_mode != 0) {
                cout << "start_mode = " << dec << start_mode
                     << " index = " << dec << index
                     << " p = " << dec << p << endl;
            }
#endif
            w256_t r2;
            if (weight_mode == max_weight_mode) {
                r2 = r;
            } else {
                for (int i = 0; i < weight_mode; i++) {
                    r2.u64[i] = r.u64[i];
                }
                for (int i = weight_mode; i < max_weight_mode; i++) {
                    r2.u64[i] = previous.u64[i];
                }
            }
#if defined(DEBUG) && 0
            cout << "w:" << dec << weight_mode
                 << " s:" << dec << start_mode << endl;
            cout << "r:" << hex << r << endl;
            cout << "p:" << hex << previous << endl;
            cout << "2:" << hex << r2 << endl;
#endif
            previous = r;
            return r2;
        }

        void setStartMode(int mode) {
            start_mode = mode;
        }

        int getStartMode() {
            return start_mode;
        }

        void setWeightMode(int mode) {
            weight_mode = mode;
            MTToolBox::setZero(previous);
        }

        int getWeightMode() {
            return weight_mode;
        }

        /**
         * @param bit_len bit length from MSB or LSB
         * @return generated numbers of bit_len
         */
        w256_t generate(int bit_len) {
            w256_t w;
#if 0
            if (reverse_bit_flag) {
                w = reverse_bit(generate());
            } else {
                w = generate();
            }
#endif
            w = generate();
            w256_t mask = make_msb_mask<w256_t>(bit_len);
            return w & mask;
        }

        /**
         * make parameters from given sequential number and
         * internal id
         * @param num sequential number
         */
        void setUpParam(ParameterGenerator& mt) {
#if defined(DEBUG)
            cout << "dSFMTAVX2 setUpParam start" << endl;
#endif
            if (size == 2) {
                param.pos1 = 1;
            } else {
                param.pos1 = mt.getUint64() % (size - 2) + 1;
            }
            if (fixedSL1 > 0) {
                param.sl1 = fixedSL1;
            } else {
                param.sl1 = mt.getUint64() % (52 - 1) + 1;
            }
            if (fixedPerm > 0) {
                param.perm = fixedPerm;
            } else {
                param.perm = (mt.getUint64() % 4) * 2 + 1;
            }
            for (int i = 0; i < 4; i++) {
                param.msk1.u64[i] = mt.getUint64() | mt.getUint64();
                param.msk1.u64[i] &= UINT64_C(0x000fffffffffffff);
            }
#if defined(DEBUG)
            cout << "dSFMTAVX2 setUpParam end" << endl;
#endif
        }

        void setZero() {
            for (int i = 0; i < size; i++) {
                MTToolBox::setZero(state[i]);
            }
            MTToolBox::setZero(lung);
            index = 0;
            MTToolBox::setZero(previous);
        }

        /**
         * This method is called by the functions in the file
         * simple_shortest_basis.hpp
         * @return true if all elements of state is zero
         */
        bool isZero() const {
            for (int i = 0; i < size; i++) {
                if (!MTToolBox::isZero(state[i])) {
                    return false;
                }
            }
            if (!MTToolBox::isZero(lung)) {
                return false;
            }
            if (weight_mode == max_weight_mode) {
                return true;
            } else {
                for (int i = max_weight_mode - weight_mode;
                     i < max_weight_mode;
                     i++) {
                    if (previous.u64[i] != 0) {
                        return false;
                    }
                }
                return true;
            }
        }

        void setParityValue(w256_t parity) {
            lung = parity;
            param.parity1 = parity;
        }

        w256_t getParityValue() const {
            return lung;
        }

        void setOneBit(int bitPos) {
            setZero();
            if (bitPos < size * element_size) {
                int idx = bitPos / element_size;
                int p = (bitPos / 52) % 4;
                int r = bitPos % 52;
                state[idx].u64[p] = UINT64_C(1) << r;
            } else {
                bitPos = bitPos - size * element_size;
                int p = (bitPos / 64) % 4;
                int r = bitPos % 64;
                lung.u64[p] = UINT64_C(1) << r;
            }
        }

        /**
         * @param that DSFMTAVX2 generator added to this generator
         */
        void add(EquidistributionCalculatable<w256_t>& other) {
            dSFMTAVX2 *that = dynamic_cast<dSFMTAVX2 *>(&other);
            if (that == 0) {
                throw std::invalid_argument(
                    "the adder should have same type as the addee.");
            }
            this->add(that);
        }

        void add(const dSFMTAVX2 * that) {
            for (int i = 0; i < size; i++) {
                state[(i + index) % size]
                    ^= that->state[(i + that->index) % size];
            }
            lung ^= that->lung;
            previous ^= that->previous;
        }

        int getMexp() const {
            return param.mexp;
        }

        int bitSize() const {
            return size * element_size + lung_size;
        }

        const std::string getHeaderString() {
            return param.get_header();
        }

        const std::string getParamString() {
            return param.get_string();
        }

        void set_reverse_bit() {
            //reverse_bit_flag = true;
        }

        void reset_reverse_bit() {
            //reverse_bit_flag = false;
        }

        void setPrefix(uint64_t p) {
            prefix = p;
        }

        void setConst() {
            const uint64_t high = UINT64_C(0x3ff0000000000000);
            dSFMTAVX2 tmp(*this);
            tmp.setZero();
            tmp.setPrefix(high);
            tmp.setup_prefix();
            tmp.next_state();
            setZero();
            setPrefix(high);
            setup_prefix();
            add(tmp);
            setPrefix(0);
        }

        bool equals(const dSFMTAVX2& that) {
            if (lung != that.lung) {
                return false;
            }
            for (int i = 0; i < size; i++) {
                if (state[(index + i) % size]
                    != that.state[(that.index + i) % size]) {
                    return false;
                }
            }
            return true;
        }

        void setFixPoint(const w256_t& fix) {
            param.fix1 = fix;
        }

        void d_p() {
            cout << "index = " << dec << index << endl;
            for (int i = 0; i < size; i++) {
                cout << state[i] << endl;
            }
            cout << lung << endl;
        }

        int periodCertification(bool noFix = false) {
            w256_t tmp;
            w256_t parity;
            parity = param.parity1;
            w256_t fix;
            fix = param.fix1;
            if (noFix) {
                tmp = lung & parity;
            } else {
                tmp = lung;
                tmp ^= fix;
                tmp = tmp & parity;
            }
            int c = count_bit(tmp) & 1;
            if (c == 1) {
                return 1;
            }
            if ((parity.u64[0] & 1) == 1) {
                lung.u64[0] ^= 1;
                return 0;
            }
            for (int i = 3; i >= 0; i--) {
                uint64_t work = 1;
                for (int j = 0; j < 64; j++) {
                    if ((work & parity.u64[i]) != 0) {
                        lung.u64[i] ^= work;
                        return 0;
                    }
                    work = work << 1;
                }
            }
            return 0;
        }
        void setFixedSL1(int value) {
            fixedSL1 = value;
        }
        void setFixedPerm(int value) {
            fixedPerm = value;
        }
    private:
        dSFMTAVX2& operator=(const dSFMTAVX2&) {
            throw std::logic_error("can't assign");
        }
        void setup_prefix() {
            const uint64_t clear = UINT64_C(0x000fffffffffffff);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 4; j++) {
                    state[i].u64[j] &= clear;
                }
            }
            if (prefix != 0) {
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < 4; j++) {
                        state[i].u64[j] |= prefix;
                    }
                }
            }
        }
        enum {sr1 = 12, lung_size = 256, element_size = 208,
              max_weight_mode = 4};
        int fixedSL1;
        int fixedPerm;
        int size;
        int index;
        int start_mode;
        int weight_mode;
        w256_t * state;
        w256_t lung;
        dSFMTAVX2_param param;
        w256_t previous;
        uint64_t prefix;
    };
}

#endif
