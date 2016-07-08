#ifndef W512_HPP
#define W512_HPP

/**
 * @file w512.hpp
 *
 * @brief 512-bit LITTLE ENDIAN class
 */

#include <iostream>
#include <iomanip>
#include <MTToolBox/util.hpp>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>

namespace MTToolBox {
    using namespace NTL;
    using namespace std;

    union w512_t {
        uint32_t u[16];
        uint64_t u64[8];
    };

    // original is 15
    static inline w512_t permutexvar_epi32(w512_t x, int perm)
    {
        w512_t r;
        for (int i = 0; i < 16; i++) {
            r.u[i] = x.u[(i + perm) % 16];
        }
        return r;
    }

    static inline w512_t SR64(w512_t x, int s)
    {
        w512_t w;
        for (int i = 0; i < 8; i++) {
            w.u64[i] = x.u64[i] >> s;
        }
        return w;
    }

    static inline w512_t SL64(w512_t x, int s)
    {
        w512_t w;
        for (int i = 0; i < 8; i++) {
            w.u64[i] = x.u64[i] << s;
        }
        return w;
    }
    /**
     * n = 1 to 512
     * upper n bit is on
     *
     */
    template<>
    inline w512_t make_msb_mask(int n)
    {
        w512_t w;
        uint64_t mask = ~UINT64_C(0);
        int n2 = 512 - n;
        int p = n2 / 64;
        for (int i = 7; i > p; i--) {
            w.u64[i] = mask;
        }
        int r = n2 % 64;
        w.u64[p] = mask << r;
        for (int i = p - 1; i >= 0; i--) {
            w.u64[i] = 0;
        }
        return w;
    }

#if 0
    static inline w512_t and_mask(w512_t a, w512_t b) {
        w512_t w;
        for (int i = 0; i < 8; i++) {
            w.u64[i] = a.u64[i] & b.u64[i];
        }
        return w;
    }
#endif

    template<>
    inline w512_t getOne() {
        w512_t one;
        one.u64[0] = 1;
        for (int i = 1; i < 8; i++) {
            one.u64[i] = 0;
        }
        return one;
    }

    /**
     * pos = 0 to 511
     *
     */
    template<>
    inline unsigned int getBitOfPos(w512_t bits, int pos) {
        if (pos >= 512) {
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
    inline void setBitOfPos(w512_t * bits, int pos, unsigned int b) {
        if (pos >= 512) {
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
    inline bool isZero(w512_t x) {
        for (int i = 0; i < 8; i++) {
            if (x.u64[i] != 0) {
                return false;
            }
        }
        return true;
    }

    template<>
    inline void setZero(w512_t& x) {
        for (int i = 0; i < 8; i++) {
            x.u64[i] = 0;
        }
    }

    inline int count_bit(w512_t x)
    {
        int c = 0;
        for (int i = 0; i < 8; i++) {
            c += count_bit(x.u64[i]);
        }
        return c;
    }

#if 0
    inline const w512_t operator>>(w512_t x, int s) {
        w512_t r;
        if (s == 0) {
            r = x;
            return r;
        } else if (s == 64) {
            r.u64[0] = x.u64[1];
            r.u64[1] = 0;
        } else if (s >= 64) {
            s = s - 64;
            r.u64[0] = x.u64[1] >> s;
            r.u64[1] = 0;
        } else {
            uint64_t w = x.u64[1] << (64 - s);
            r.u64[1] = x.u64[1] >> s;
            r.u64[0] = w | (x.u64[0] >> s);
        }
        return r;
    }

    inline const w512_t operator<<(w512_t x, int s) {
        w512_t r;
        if (s == 0) {
            r = x;
            return r;
        } else if (s == 64) {
            r.u64[1] = x.u64[0];
            r.u64[0] = 0;
        } else if (s > 64) {
            s = s - 64;
            r.u64[1] = x.u64[0] << s;
            r.u64[0] = 0;
        } else {
            uint64_t w = x.u64[0] >> (64 - s);
            r.u64[0] = x.u64[0] << s;
            r.u64[1] = w | (x.u64[1] << s);
        }
        return r;
    }
#endif

    inline const w512_t operator&(w512_t x, w512_t y) {
        w512_t r;
        for (int i = 0; i < 8; i++) {
            r.u64[i] = x.u64[i] & y.u64[i];
        }
        return r;
    }

#if 0
    inline const w512_t operator^(w512_t x, w512_t y) {
        w512_t r = x;
        r.u64[0] ^= y.u64[0];
        r.u64[1] ^= y.u64[1];
        return r;
    }

    inline const w512_t operator~(w512_t x) {
        w512_t r;
        for (int i = 0; i < 8; i++) {
            r.u64[i] = ~x.u64[i];
        }
        return r;
    }
#endif

    inline w512_t& operator|=(w512_t& x, w512_t y) {
        for (int i = 0; i < 8; i++) {
            x.u64[i] |= y.u64[i];
        }
        return x;
    }

    inline w512_t& operator^=(w512_t& x, w512_t y) {
        for (int i = 0; i < 8; i++) {
            x.u64[i] ^= y.u64[i];
        }
        return x;
    }

    inline bool operator==(const w512_t& x, const w512_t y) {
        for (int i = 0; i < 8; i++) {
            if (x.u64[i] != y.u64[i]) {
                return false;
            }
        }
        return true;
    }

    inline bool operator!=(const w512_t& x, const w512_t y) {
        for (int i = 0; i < 8; i++) {
            if (x.u64[i] != y.u64[i]) {
                return true;
            }
        }
        return false;
    }

    inline ostream& operator<<(ostream& os, w512_t x) {
        for (int i = 7; i >= 0; i--) {
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

    inline istream& operator>>(istream&is, w512_t& x) {
        setZero(x);
        char buff[150];
        int length = 136;
        buff[136] = 0;
        is.read(buff, length);
        char * p = buff;
        for (int i = 7; i >= 0; i--) {
            x.u64[i] = strtoull(p, NULL, 16);
            p += 17;
        }
        return is;
    }

    static inline int calc_1pos(w512_t x)
    {
        if (isZero(x)) {
            return -1;
        }
        int p = 0;
        for (int i = 0; i < 8; i++) {
            if (x.u64[i] != 0) {
                p = i;
                break;
            }
        }
        int64_t y = (int64_t)x.u64[p];
        y = count_bit((uint64_t)(y & -y) - 1);
        return 511 - y - p * 64;
    }

    template<>
    inline w512_t convert(uint64_t x) {
        w512_t w;
        w.u64[0] = x;
        for (int i = 1; i < 8; i++) {
            w.u64[i] = 0;
        }
        return w;
    }

    template<>
    inline w512_t convert(uint32_t x) {
        w512_t w;
        w.u64[0] = x;
        for (int i = 1; i < 8; i++) {
            w.u64[i] = 0;
        }
        return w;
    }

    static inline w512_t reverse_bit(w512_t x) {
        w512_t w;
        for (int i = 0; i < 8; i++) {
            w.u64[i] = reverse_bit(x.u64[7 - i]);
        }
        return w;
    }
}
#endif
