#pragma once
#ifndef W256_HPP
#define W256_HPP
/**
 * @file w256.hpp
 *
 * @brief 256-bit LITTLE ENDIAN class
 *
 */

#include "devavxprng.h"

namespace MTToolBox {
    using namespace NTL;
    using namespace std;

    union w256_t {
        uint32_t u[8];
        uint64_t u64[4];
    };

    static inline w256_t permutexvar_epi32(w256_t x, int perm)
    {
        w256_t r;
        for (int i = 0; i < 8; i++) {
            r.u[i] = x.u[(i + perm) % 8];
        }
        return r;
    }

    static inline w256_t SR64(w256_t x, int s)
    {
        w256_t w;
        for (int i = 0; i < 4; i++) {
            w.u64[i] = x.u64[i] >> s;
        }
        return w;
    }

    static inline w256_t SL64(w256_t x, int s)
    {
        w256_t w;
        for (int i = 0; i < 4; i++) {
            w.u64[i] = x.u64[i] << s;
        }
        return w;
    }

    /**
     * n = 1 to 256
     * upper n bit is on
     *
     */
    static inline w256_t make_msb_mask(int n)
    {
        w256_t w;
        uint64_t mask = ~UINT64_C(0);
        int n2 = 256 - n;
        int p = n2 / 64;
        for (int i = 3; i > p; i--) {
            w.u64[i] = mask;
        }
        int r = n2 % 64;
        w.u64[p] = mask << r;
        for (int i = p - 1; i >= 0; i--) {
            w.u64[i] = 0;
        }
        return w;
    }

    template<>
    inline w256_t getOne() {
        w256_t one;
        one.u64[0] = 1;
        for (int i = 1; i < 4; i++) {
            one.u64[i] = 0;
        }
        return one;
    }

    /**
     * pos = 0 to 511
     *
     */
    template<>
    inline unsigned int getBitOfPos(w256_t bits, int pos) {
        if (pos >= 256) {
            return 0;
        }
        int p = pos / 64;
        int r = pos % 64;
        return (bits.u64[p] >> r) & 1;
    }

    /**
     * pos = 0 to 511
     * if b is 0, bit will be reset.
     */
    template<>
    inline void setBitOfPos(w256_t * bits, int pos, unsigned int b) {
        if (pos >= 256) {
            return;
        }
        int p = pos / 64;
        int r = pos % 64;
        uint64_t x = b & 1;
        uint64_t mask = ~UINT64_C(1) << r;
        bits->u64[p] &= mask;
        bits->u64[p] |= x << r;
    }

    template<>
    inline bool isZero(w256_t x) {
        for (int i = 0; i < 4; i++) {
            if (x.u64[i] != 0) {
                return false;
            }
        }
        return true;
    }

    template<>
    inline void setZero(w256_t& x) {
        for (int i = 0; i < 4; i++) {
            x.u64[i] = 0;
        }
    }

    inline int count_bit(w256_t x)
    {
        int c = 0;
        for (int i = 0; i < 4; i++) {
            c += count_bit(x.u64[i]);
        }
        return c;
    }

    inline const w256_t operator&(w256_t x, w256_t y)
    {
        w256_t r;
        for (int i = 0; i < 4; i++) {
            r.u64[i] = x.u64[i] & y.u64[i];
        }
        return r;
    }

    inline w256_t& operator|=(w256_t& x, w256_t y)
    {
        for (int i = 0; i < 4; i++) {
            x.u64[i] |= y.u64[i];
        }
        return x;
    }

    inline w256_t& operator^=(w256_t& x, w256_t y)
    {
        for (int i = 0; i < 4; i++) {
            x.u64[i] ^= y.u64[i];
        }
        return x;
    }

    inline bool operator==(const w256_t& x, const w256_t y)
    {
        for (int i = 0; i < 4; i++) {
            if (x.u64[i] != y.u64[i]) {
                return false;
            }
        }
        return true;
    }

    inline bool operator!=(const w256_t& x, const w256_t y)
    {
        for (int i = 0; i < 4; i++) {
            if (x.u64[i] != y.u64[i]) {
                return true;
            }
        }
        return false;
    }

    inline ostream& operator<<(ostream& os, w256_t x)
    {
        for (int i = 3; i >= 0; i--) {
            os << hex;
            os << setfill('0');
            os << setw(16);
            os << x.u64[i];
            if (i != 0) {
                os << ".";
            }
        }
        return os;
    }

    static inline int calc_1pos(w256_t x)
    {
        if (isZero(x)) {
            return -1;
        }
        int p = 0;
        for (int i = 0; i < 4; i++) {
            if (x.u64[i] != 0) {
                p = i;
                break;
            }
        }
        int64_t y = (int64_t)x.u64[p];
        y = count_bit((uint64_t)(y & -y) - 1);
        return 255 - y - p * 64;
    }

    template<>
    inline w256_t convert(uint64_t x)
    {
        w256_t w;
        w.u64[0] = x;
        for (int i = 1; i < 4; i++) {
            w.u64[i] = 0;
        }
        return w;
    }

    template<>
    inline w256_t convert(uint32_t x)
    {
        w256_t w;
        w.u64[0] = x;
        for (int i = 1; i < 4; i++) {
            w.u64[i] = 0;
        }
        return w;
    }

    static inline w256_t reverse_bit(w256_t x)
    {
        w256_t w;
        for (int i = 0; i < 4; i++) {
            w.u64[i] = reverse_bit(x.u64[3 - i]);
        }
        return w;
    }
}
#endif
