#pragma once
#ifndef SFMTAVX512FSEARCH_HPP
#define SFMTAVX512FSEARCH_HPP
/**
 * @file SFMTAVX512Fsearch.hpp
 */

#include "devavxprng.h"
#include "w512.hpp"
#include <sstream>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    using namespace NTL;
    using namespace std;

    class SFMTAVX512F_param {
    public:
        int mexp;
        int pos1;
        int sl1;
        int sr1;
        int perm;
        w512_t mat1;
        w512_t parity1;

        SFMTAVX512F_param() {
            mexp = 0;
            pos1 = 0;
            sl1 = 0;
            sr1 = 0;
            perm = 1;
            MTToolBox::setZero(mat1);
            MTToolBox::setZero(parity1);
        }

        SFMTAVX512F_param(const SFMTAVX512F_param& src) {
            mexp = src.mexp;
            pos1 = src.pos1;
            sl1 = src.sl1;
            sr1 = src.sr1;
            perm = src.perm;
            mat1 = src.mat1;
            parity1 = src.parity1;
        }

        /**
         * This method is used in output.hpp.
         * @return header line of output.
         */
        const string get_header() const {
            return "mexp, pos1, sl1, sr1, perm, mat1, parity1";
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
            ss << dec << sr1 << ",";
            ss << dec << perm << ",";
            ss << mat1 << ",";
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
            ss << "sr1:" << dec << sr1 << endl;
            ss << "perm:" << dec << perm << endl;
            ss << "mat1:" << mat1 << endl;
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
            sr1 = strtoul(para, &para, 10);
            para++;
            perm = strtoul(para, &para, 10);
            para++;
            for (int i = 7; i >= 0; i--) {
                mat1.u64[i] = strtoull(para, &para, 16);
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
     * @class SFMTAVX512F
     * @brief SFMTAVX512F generator class used for dynamic creation
     *
     * This class is one of the main class of SFMTAVX512F dynamic creator.
     * This class is designed to be called from programs in MTToolBox,
     * but is not a subclass of some abstract class.
     * Instead, this class is passed to them as template parameters.
     */
    class SFMTAVX512F : public ReducibleGenerator<w512_t, uint32_t> {
    public:
        /**
         * Constructor by mexp.
         * @param mexp Mersenne Exponent
         */
        SFMTAVX512F(int mexp) {
            size = mexp / 512;
            state = new w512_t[size];
            param.mexp = mexp;
            param.pos1 = 0;
            param.sl1 = 0;
            param.sr1 = 0;
            param.perm = 1;
            fixed = false;
            MTToolBox::setZero(param.mat1);
            MTToolBox::setZero(param.parity1);
            MTToolBox::setZero(lung);
            index = 0;
            reverse_bit_flag = false;
            start_mode = 0;
            weight_mode = 16;
            MTToolBox::setZero(previous);
        }

        ~SFMTAVX512F() {
            delete[] state;
        }

        /**
         * The copy constructor.
         * @param src The origin of copy.
         */
        SFMTAVX512F(const SFMTAVX512F& src) : param(src.param) {
            size = src.size;
            state = new w512_t[size];
            for (int i = 0; i < size; i++) {
                state[i] = src.state[i];
            }
            fixed = src.fixed;
            lung = src.lung;
            index = src.index;
            start_mode = src.start_mode;
            weight_mode = src.weight_mode;
            reverse_bit_flag = src.reverse_bit_flag;
            previous = src.previous;
        }

        /**
         * Constructor by parameter.
         * @param src_param
         */
        SFMTAVX512F(const SFMTAVX512F_param& src_param) : param(src_param) {
            size = src_param.mexp / 512;
            state = new w512_t[size];
            for (int i = 0; i < size; i++) {
                MTToolBox::setZero(state[i]);
            }
            MTToolBox::setZero(lung);
            fixed = false;
            index = 0;
            start_mode = 0;
            weight_mode = 16;
            MTToolBox::setZero(previous);
            reverse_bit_flag = false;
        }

        EquidistributionCalculatable<w512_t, uint32_t> * clone() const {
            return new SFMTAVX512F(*this);
        }

        /**
         * This method initialize internal state.
         * This initialization is simple.
         * @param seed seed for initialization
         */
        void seed(w512_t seed) {
            setZero();
            state[0] = seed;
            uint32_t * pstate = new uint32_t[(size + 1) * 16];
            for (int i = 0; i < 16; i++) {
                pstate[i] = seed.u[i];
            }
            for (int i = 16; i < (size + 1) * 16; i++) {
                pstate[i] = 0;
            }
            for (int i = 1; i < size * 16; i++) {
                pstate[i] ^= i + UINT32_C(1812433253)
                    * (pstate[i - 1] ^ (pstate[i - 1] >> 30));
            }
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 16; j++) {
                    state[i].u[j] = pstate[i * 16 + j];
                }
            }
            for (int i = 0; i < 16; i++) {
                lung.u[i] = pstate[size * 16 + i];
            }
            index = 0;
            delete[] pstate;
            MTToolBox::setZero(previous);
        }

        void do_recursion(w512_t *r, w512_t *a, w512_t *b, w512_t *lung) {
            *lung = permutexvar_epi32(*lung, param.perm);
            *lung ^= *b;
            *lung ^= SL64(*a, param.sl1);
            w512_t t;
            t = *lung;
            t ^= SR64(*b & param.mat1, param.sr1);
            t ^= *a;
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
                for (int i = 0; i < 16 - start_mode; i++) {
                    r.u[i] = state[p].u[i + start_mode];
                }
                int j = 0;
                for (int i = 16 - start_mode; i < 16; i++) {
                    r.u[i] = state[index].u[j++];
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
                    r2.u[i] = r.u[i];
                }
                for (int i = weight_mode; i < max_weight_mode; i++) {
                    r2.u[i] = previous.u[i];
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
            if (reverse_bit_flag) {
                w = reverse_bit(generate());
            } else {
                w = generate();
            }
            w512_t mask = make_msb_mask<w512_t>(bit_len);
            return w & mask;
        }

        /**
         * make parameters from given sequential number and
         * internal id
         * @param num sequential number
         */
        void setUpParam(AbstractGenerator<uint32_t>& mt) {
            if (param.mexp == 1279) {
                param.pos1 = 1;
            } else {
                param.pos1 = mt.generate() % (size - 2) + 1;
            }
            if (fixed) {
                // These parameters are not best ones.
                param.sl1 = 19;
                param.sr1 = 7;
            } else {
                param.sl1 = mt.generate() % (64 - 1) + 1;
                param.sr1 = mt.generate() % (64 - 1) + 1;
            }
            param.perm = (mt.generate() % 8) * 2 + 1;
            for (int i = 0; i < 16; i++) {
                param.mat1.u[i] = mt.generate() | mt.generate();
            }
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
                    if (previous.u[i] != 0) {
                        return false;
                    }
                }
            }
            return true;
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
                int p = (bitPos / 64) % 8;
                int r = bitPos % 64;
                state[idx].u64[p] = UINT64_C(1) << r;
            } else {
                bitPos = bitPos - size * element_size;
                int p = (bitPos / 64) % 8;
                int r = bitPos % 64;
                lung.u64[p] = UINT64_C(1) << r;
            }
        }

        /**
         * This method is called by functions in the file
         * simple_shortest_basis.hpp addition of internal state as
         * GF(2) vector is possible when state transition function and
         * output function is GF(2)-linear.
         * @param that SFMTAVX512F generator added to this generator
         */
        void add(EquidistributionCalculatable<w512_t, uint32_t>& other) {
            SFMTAVX512F *that = dynamic_cast<SFMTAVX512F *>(&other);
            if (that == 0) {
                throw std::invalid_argument(
                    "the adder should have same type as the addee.");
            }
            this->add(that);
        }

        void add(const SFMTAVX512F * that) {
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
            return (size + 1) * 512;
        }

        const std::string getHeaderString() {
            return param.get_header();
        }

        const std::string getParamString() {
            return param.get_string();
        }

        /**
         * This method is called by the functions in search_temper.hpp
         * to calculate the equidistribution properties from LSB
         */
        void set_reverse_bit() {
            reverse_bit_flag = true;
        }

        /**
         * This method is called by the functions in search_temper.hpp
         * to reset the reverse_bit_flag
         */
        void reset_reverse_bit() {
            reverse_bit_flag = false;
        }
        int periodCertification() {
            w512_t tmp = lung & param.parity1;
            int c = count_bit(tmp) & 1;
            if ((c & 1) == 1) {
                return 1;
            }
            if ((param.parity1.u64[0] & 1) == 1) {
                lung.u64[0] ^= 1;
                return 0;
            }
            for (int i = 7; i >= 0; i--) {
                uint64_t work = 1;
                for (int j = 0; j < 64; j++) {
                    if ((work & param.parity1.u[i]) != 0) {
                        lung.u64[i] ^= work;
                        return 0;
                    }
                    work = work << 1;
                }
            }
            return 0;
        }

        void d_p() {
            cout << "index = " << dec << index << endl;
            for (int i = 0; i < size; i++) {
                cout << state[i] << endl;
            }
            cout << lung << endl;
        }
        void setFixed(bool value) {
            fixed = value;
        }

    private:
        SFMTAVX512F& operator=(const SFMTAVX512F&) {
            throw std::logic_error("can't assign");
        }
        enum {element_size = 512, max_weight_mode = 16};
        bool fixed;
        int size;
        int index;
        int start_mode;
        int weight_mode;
        w512_t * state;
        w512_t lung;
        SFMTAVX512F_param param;
        bool reverse_bit_flag;
        w512_t previous;
    };
}

#endif
