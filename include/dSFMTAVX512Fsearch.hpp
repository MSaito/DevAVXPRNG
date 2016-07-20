#pragma once
#ifndef DSFMTAVX512FSEARCH_HPP
#define DSFMTAVX512FSEARCH_HPP
/**
 * @file dSFMTAVX512Fsearch.hpp
 *
 */

#include "devavxprng.h"
#include <sstream>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/MersenneTwister64.hpp>
#include <MTToolBox/util.hpp>
#include "w512.hpp"

/**
 * @namespace dSFMTAVX512F
 * name space for dSFMTAVX512F
 */
namespace MTToolBox {
    using namespace NTL;
    using namespace std;

    /**
     * @class dSFMTAVX512F_param
     * @brief a class keeping parameters of dSFMTAVX512F
     *
     * This class keeps parameters of dSFMTAVX512F, and has some
     * method for outputting parameters.
     */
    class dSFMTAVX512F_param {
    public:
        int mexp;
        int pos1;
        int sl1;
        int sl2;
        int perm;
        w512_t msk1;
        w512_t fix1;
        w512_t parity1;

        dSFMTAVX512F_param() {
            mexp = 0;
            pos1 = 0;
            sl1 = 0;
            sl2 = 0;
            perm = 1;
            setZero(msk1);
            setZero(fix1);
            setZero(parity1);
        }

        dSFMTAVX512F_param(const dSFMTAVX512F_param& src) {
            mexp = src.mexp;
            pos1 = src.pos1;
            sl1 = src.sl1;
            sl2 = src.sl2;
            perm = src.perm;
            msk1 = src.msk1;
            fix1 = src.fix1;
            parity1 = src.parity1;
        }
        /**
         * This method is used in output.hpp.
         * @return header line of output.
         */
        const string get_header() const {
            return "mexp, pos1, sl1, sl2, perm, msk, fix, parity";
        }

        /**
         * This method is used in output.hpp.
         * @return string of parameters
         */
        const string get_string() const {
            stringstream ss;
            ss << dec << mexp << ",";
            ss << dec << pos1 << ",";
            ss << dec << sl1 << ",";
            ss << dec << sl2 << ",";
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
            ss << "sl2:" << dec << sl2 << endl;
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
            sl2 = strtoul(para, &para, 10);
            para++;
            perm = strtoul(para, &para, 10);
            para++;
            for (int i = 7; i >= 0; i--) {
                msk1.u64[i] = strtoull(para, &para, 16);
                para++;
            }
            if (para >= str + len) {
                return;
            }
            for (int i = 7; i >= 0; i--) {
                fix1.u64[i] = strtoull(para, &para, 16);
                para++;
            }
            if (para >= str + len) {
                return;
            }
            for (int i = 7; i >= 0; i--) {
                parity1.u64[i] = strtoull(para, &para, 16);
                para++;
            }
        }
    };

    /**
     * @class dSFMTAVX512F
     * @brief DSFMTAVX512F generator class used for dynamic creation
     *
     * This class is one of the main class of DSFMTAVX512F dynamic creator.
     * This class is designed to be called from programs in MTToolBox,
     * but is not a subclass of some abstract class.
     * Instead, this class is passed to them as template parameters.
     */
    class dSFMTAVX512F : public ReducibleGenerator<w512_t, uint64_t> {
    public:
        /**
         * Constructor by mexp.
         * @param mexp Mersenne Exponent
         */
        dSFMTAVX512F(int mexp) {
#if defined(DEBUG)
            cout << "dSFMTAVX512F constructor start" << endl;
#endif
            size = (mexp - lung_size) / element_size + 1;
            state = new w512_t[size];
            param.mexp = mexp;
            param.pos1 = 0;
            param.sl1 = 0;
            param.sl2 = 0;
            param.perm = 1;
            MTToolBox::setZero(param.msk1);
            MTToolBox::setZero(param.parity1);
            index = 0;
            start_mode = 0;
            weight_mode = max_weight_mode;
            MTToolBox::setZero(previous);
            MTToolBox::setZero(lung);
            prefix = 0;
            fixed = false;
            fixedSL1 = 0;
        }

        ~dSFMTAVX512F() {
            delete[] state;
        }

        /**
         * The copy constructor.
         * @param src The origin of copy.
         */
        dSFMTAVX512F(const dSFMTAVX512F& src) : param(src.param) {
#if defined(DEBUG)
            cout << "dSFMTAVX512F constructor start" << endl;
#endif
            size = src.size;
            state = new w512_t[size];
            for (int i = 0; i < size; i++) {
                state[i] = src.state[i];
            }
            lung = src.lung;
            index = src.index;
            start_mode = src.start_mode;
            weight_mode = src.weight_mode;
            previous = src.previous;
            prefix = src.prefix;
            fixed = src.fixed;
            fixedSL1 = src.fixedSL1;
        }

        /**
         * Constructor by parameter.
         * @param src_param
         */
        dSFMTAVX512F(const dSFMTAVX512F_param& src_param) : param(src_param) {
#if defined(DEBUG)
            cout << "dSFMTAVX512F constructor start" << endl;
#endif
            size = (src_param.mexp - lung_size) / element_size + 1;
            state = new w512_t[size];
            for (int i = 0; i < size; i++) {
                MTToolBox::setZero(state[i]);
            }
            index = 0;
            start_mode = 0;
            weight_mode = max_weight_mode;
            MTToolBox::setZero(previous);
            MTToolBox::setZero(lung);
            prefix = 0;
            fixed = false;
            fixedSL1 = 0;
        }

        EquidistributionCalculatable<w512_t, uint64_t> * clone() const {
            return new dSFMTAVX512F(*this);
        }

        /**
         * This method initialize internal state.
         * This initialization is simple.
         * @param seed seed for initialization
         */
        void seed(w512_t seed) {
            //setZero();
            uint64_t * pstate = new uint64_t[(size + 1) * 8];
            for (int i = 0; i < (size + 1) * 8; i++) {
                pstate[i] = 0;
            }
            for (int i = 0; i < 8; i++) {
                pstate[i] = seed.u64[i];
            }
            for (int i = 1; i < (size + 1) * 8; i++) {
                pstate[i] ^= i + UINT64_C(6364136223846793005)
                    * (pstate[i - 1] ^ (pstate[i - 1] >> 62));
            }
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 8; j++) {
                    state[i].u64[j] = pstate[i * 8 + j];
                }
            }
            for (int i = 0; i < 8; i++) {
                lung.u64[i] = pstate[size * 8 + i];
            }
            delete[] pstate;
            index = 0;
            MTToolBox::setZero(previous);
            setup_prefix();
        }

        void do_recursion(w512_t *r, w512_t *a, w512_t *b, w512_t *lung) {
            *lung = permutexvar_epi32(*lung, param.perm);
            *lung ^= *b;
            *lung ^= SL64(*a, param.sl1);
            w512_t t;
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
        w512_t generate() {
            next_state();
            w512_t r;
            index = index % size;
            int p = (index + size - 1) % size;
            if (start_mode == 0) {
                r = state[index];
            } else {
                for (int i = 0; i < 8 - start_mode; i++) {
                    r.u64[i] = state[p].u64[i + start_mode];
                }
                int j = 0;
                for (int i = 8 - start_mode; i < 8; i++) {
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
            w512_t r2;
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

        w512_t generate(int bit_len) {
            w512_t w;
            w = generate();
            w512_t mask = make_msb_mask<w512_t>(bit_len);
            return w & mask;
        }

        void setUpParam(AbstractGenerator<uint64_t>& mt) {
#if defined(DEBUG)
            cout << "dSFMTAVX512Fsearch setUpParam start" << endl;
#endif
            if (size == 2) {
                param.pos1 = 1;
            } else {
                param.pos1 = mt.generate() % (size - 2) + 1;
            }
            if (fixed) {
                param.sl1 = fixedSL1;
                param.sl2 = 0;
            } else {
                param.sl1 = mt.generate() % (52 - 1) + 1;
                param.sl2 = 0;
            }
            param.perm = (mt.generate() % 8) * 2 + 1;
            for (int i = 0; i < 8; i++) {
                param.msk1.u64[i] = mt.generate() | mt.generate();
                param.msk1.u64[i] &= UINT64_C(0x000fffffffffffff);
            }
#if defined(DEBUG)
            cout << "dSFMTAVX512F setUpParam end" << endl;
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
                for (int i = max_weight_mode - weight_mode; i < max_weight_mode;
                     i++) {
                    if (previous.u64[i] != 0) {
                        return false;
                    }
                }
                return true;
            }
        }

        void setParityValue(w512_t parity) {
            lung = parity;
            param.parity1 = parity;
        }

        w512_t getParityValue() const {
            return lung;
        }

        void setOneBit(int bitPos) {
            setZero();
            if (bitPos < size * element_size) {
                int idx = bitPos / element_size;
                int p = (bitPos / 52) % 8;
                int r = bitPos % 52;
                state[idx].u64[p] = UINT64_C(1) << r;
            } else {
                bitPos = bitPos - size * element_size;
                int p = (bitPos / 64) % 8;
                int r = bitPos % 64;
                lung.u64[p] = UINT64_C(1) << r;
            }
        }

        void add(EquidistributionCalculatable<w512_t, uint64_t>& other) {
            dSFMTAVX512F *that = dynamic_cast<dSFMTAVX512F *>(&other);
            if (that == 0) {
                throw std::invalid_argument(
                    "the adder should have same type as the addee.");
            }
            this->add(that);
        }

        void add(const dSFMTAVX512F * that) {
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
            dSFMTAVX512F tmp(*this);
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

        bool equals(const dSFMTAVX512F& that) {
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

        void setFixPoint(const w512_t& fix) {
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
            w512_t tmp;
            w512_t parity;
            parity = param.parity1;
            w512_t fix;
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
            for (int i = 7; i >= 0; i--) {
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
        void setFixed(bool value) {
            fixed = value;
        }
        void setFixedSL1(int value) {
            fixedSL1 = value;
        }
    private:
        dSFMTAVX512F& operator=(const dSFMTAVX512F&) {
            throw std::logic_error("can't assign");
        }
        void setup_prefix() {
            const uint64_t clear = UINT64_C(0x000fffffffffffff);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 8; j++) {
                    state[i].u64[j] &= clear;
                }
            }
            if (prefix != 0) {
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < 8; j++) {
                        state[i].u64[j] |= prefix;
                    }
                }
            }
        }
        enum {sr1 = 12, lung_size = 512, element_size = 416,
              max_weight_mode = 8};
        int fixedSL1;
        bool fixed;
        int size;
        int index;
        int start_mode;
        int weight_mode;
        w512_t * state;
        w512_t lung;
        dSFMTAVX512F_param param;
        w512_t previous;
        uint64_t prefix;
    };
}

#endif
