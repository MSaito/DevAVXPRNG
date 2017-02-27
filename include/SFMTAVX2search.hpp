#pragma once
#ifndef SFMTAVX2SEARCH_HPP
#define SFMTAVX2SEARCH_HPP
/**
 * @file SFMTAVX2search.hpp
 */

#include "devavxprng.h"
#include "w256.hpp"
#include <sstream>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    using namespace NTL;
    using namespace std;

    class SFMTAVX2_param {
    public:
        int mexp;
        int pos1;
        int sl1;
        int sr1;
        int perm;
        w256_t mat1;
        w256_t parity1;

        SFMTAVX2_param() {
            mexp = 0;
            pos1 = 0;
            sl1 = 0;
            sr1 = 0;
            perm = 1;
            MTToolBox::setZero(mat1);
            MTToolBox::setZero(parity1);
        }

        SFMTAVX2_param(const SFMTAVX2_param& src) {
            mexp = src.mexp;
            pos1 = src.pos1;
            sl1 = src.sl1;
            sr1 = src.sr1;
            perm = src.perm;
            mat1 = src.mat1;
            parity1 = src.parity1;
        }

        const string get_header() const {
            return "mexp, pos1, sl1, sr1, perm, mat1, parity1";
        }

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
            for (int i = 3; i >= 0; i--) {
                mat1.u64[i] = strtoull(para, &para, 16);
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

    class SFMTAVX2 : public ReducibleGenerator<w256_t> {
    public:
        /**
         * Constructor by mexp.
         * @param mexp Mersenne Exponent
         */
        SFMTAVX2(int mexp) {
            size = mexp / 256;
            state = new w256_t[size];
            param.mexp = mexp;
            param.pos1 = 0;
            param.sl1 = 0;
            param.sr1 = 0;
            param.perm = 1;
            fixedSL1 = 0;
            fixedSR1 = 0;
            fixedPerm = 0;
            MTToolBox::setZero(param.mat1);
            MTToolBox::setZero(param.parity1);
            MTToolBox::setZero(lung);
            index = 0;
            reverse_bit_flag = false;
            start_mode = 0;
            weight_mode = max_weight_mode;
            MTToolBox::setZero(previous);
        }

        ~SFMTAVX2() {
            delete[] state;
        }

        SFMTAVX2(const SFMTAVX2& src) : param(src.param) {
            size = src.size;
            state = new w256_t[size];
            for (int i = 0; i < size; i++) {
                state[i] = src.state[i];
            }
            lung = src.lung;
            index = src.index;
            start_mode = src.start_mode;
            weight_mode = src.weight_mode;
            fixedSL1 = src.fixedSL1;
            fixedSR1 = src.fixedSR1;
            fixedPerm = src.fixedPerm;
            reverse_bit_flag = src.reverse_bit_flag;
            previous = src.previous;
        }

        SFMTAVX2(const SFMTAVX2_param& src_param) : param(src_param) {
            size = src_param.mexp / 256;
            state = new w256_t[size];
            for (int i = 0; i < size; i++) {
                MTToolBox::setZero(state[i]);
            }
            MTToolBox::setZero(lung);
            index = 0;
            start_mode = 0;
            fixedSL1 = 0;
            fixedSR1 = 0;
            fixedPerm = 0;
            weight_mode = max_weight_mode;
            MTToolBox::setZero(previous);
            reverse_bit_flag = false;
        }

        EquidistributionCalculatable<w256_t> * clone() const {
            return new SFMTAVX2(*this);
        }

        void seed(w256_t seed) {
            setZero();
            state[0] = seed;
            uint32_t * pstate = new uint32_t[(size + 1) * 8];
            for (int i = 0; i < 8; i++) {
                pstate[i] = seed.u[i];
            }
            for (int i = 8; i < (size + 1) * 8; i++) {
                pstate[i] = 0;
            }
            for (int i = 1; i < (size + 1) * 8; i++) {
                pstate[i] ^= i + UINT32_C(1812433253)
                    * (pstate[i - 1] ^ (pstate[i - 1] >> 30));
            }
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 8; j++) {
                    state[i].u[j] = pstate[i * 8 + j];
                }
            }
            for (int i = 0; i < 8; i++) {
                lung.u[i] = pstate[size * 8 + i];
            }
            index = 0;
            delete[] pstate;
            MTToolBox::setZero(previous);
        }

        void do_recursion(w256_t *r, w256_t *a, w256_t *b, w256_t *lung) {
            *lung = permutexvar_epi32(*lung, param.perm);
            *lung ^= *b;
            *lung ^= SL64(*a, param.sl1);
            w256_t t;
            t = *lung;
            t ^= SR64(*b & param.mat1, param.sr1);
            t ^= *a;
            *r = t;
        }

        void next_state() {
            index = (index + 1) % size;
            do_recursion(&state[index],
                         &state[index],
                         &state[(index + param.pos1) % size],
                         &lung);
        }

        w256_t generate() {
            next_state();
            w256_t r;
            index = index % size;
            int p = (index + size - 1) % size;
            if (start_mode == 0) {
                r = state[index];
            } else {
                for (int i = 0; i < 8 - start_mode; i++) {
                    r.u[i] = state[p].u[i + start_mode];
                }
                int j = 0;
                for (int i = 8 - start_mode; i < 8; i++) {
                    r.u[i] = state[index].u[j++];
                }
            }
            w256_t r2;
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

        w256_t generate(int bit_len) {
            w256_t w;
            if (reverse_bit_flag) {
                w = reverse_bit(generate());
            } else {
                w = generate();
            }
            w256_t mask = make_msb_mask<w256_t>(bit_len);
            return w & mask;
        }

        void setUpParam(ParameterGenerator& mt) {
            if (size == 2) {
                param.pos1 = 1;
            } else {
                param.pos1 = mt.getUint32() % (size - 1) + 1;
            }
            if (fixedSL1 > 0) {
                param.sl1 = fixedSL1;
            } else {
                param.sl1 = mt.getUint32() % 32 + 1;
            }
            if (fixedSR1 > 0) {
                param.sr1 = fixedSR1;
            } else {
                param.sr1 = mt.getUint32() % 32 + 1;
            }
            if (fixedPerm > 0) {
                param.perm = fixedPerm;
            } else {
                param.perm = (mt.getUint32() % 4) * 2 + 1;
            }
            for (int i = 0; i < 8; i++) {
                param.mat1.u[i] = mt.getUint32() | mt.getUint32();
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
                int p = (bitPos / 64) % 4;
                int r = bitPos % 64;
                state[idx].u64[p] = UINT64_C(1) << r;
            } else {
                bitPos = bitPos - size * element_size;
                int p = (bitPos / 64) % 4;
                int r = bitPos % 64;
                lung.u64[p] = UINT64_C(1) << r;
            }
        }

        void add(EquidistributionCalculatable<w256_t>& other) {
            SFMTAVX2 *that = dynamic_cast<SFMTAVX2 *>(&other);
            if (that == 0) {
                throw std::invalid_argument(
                    "the adder should have same type as the addee.");
            }
            this->add(that);
        }

        void add(const SFMTAVX2 * that) {
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
            return (size + 1) * 256;
        }

        const std::string getHeaderString() {
            return param.get_header();
        }

        const std::string getParamString() {
            return param.get_string();
        }

        void set_reverse_bit() {
            reverse_bit_flag = true;
        }

        void reset_reverse_bit() {
            reverse_bit_flag = false;
        }

        int periodCertification() {
            w256_t tmp = lung & param.parity1;
            int c = count_bit(tmp) & 1;
            if ((c & 1) == 1) {
                return 1;
            }
            if ((param.parity1.u64[0] & 1) == 1) {
                lung.u64[0] ^= 1;
                return 0;
            }
            for (int i = 3; i >= 0; i--) {
                uint64_t work = 1;
                for (int j = 0; j < 64; j++) {
                    if ((work & param.parity1.u64[i]) != 0) {
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
        void setFixedSL1(int value) {
            fixedSL1 = value;
        }
        void setFixedSR1(int value) {
            fixedSR1 = value;
        }
        void setFixedPerm(int value) {
            fixedPerm = value;
        }

    private:
        SFMTAVX2& operator=(const SFMTAVX2&) {
            throw std::logic_error("can't assign");
        }
        enum {element_size = 256, max_weight_mode = 8};
        int fixedSL1;
        int fixedSR1;
        int fixedPerm;
        int size;
        int index;
        int start_mode;
        int weight_mode;
        w256_t * state;
        w256_t lung;
        SFMTAVX2_param param;
        bool reverse_bit_flag;
        w256_t previous;
    };
}

#endif
