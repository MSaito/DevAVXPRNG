#pragma once
#ifndef MTTOOLBOX_ALGORITHM_SIMD_EQUIDISTRIBUTION_SIMPLE_HPP
#define MTTOOLBOX_ALGORITHM_SIMD_EQUIDISTRIBUTION_SIMPLE_HPP
/**
 * @file AlgorithmSIMDEquidistribution.hpp
 */

#include "devavxprng.h"

#if HAVE_STD_SP
#include <memory>
#elif HAVE_STD_TR1_SP
#include <tr1/memory>
#else
#pragma GCC error "do not have shared_ptr"
#endif

#include <stdexcept>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
#if HAVE_STD_SP
    using std::shared_ptr;
#else
    using std::tr1::shared_ptr;
#endif

    struct SIMDInfo {
        bool fastMode;
        int bitSize;
        int bitMode;
        int elementNo;
    };

    /**
     * @class linear_generator_vector
     *\japanese
     * @brief GF(2)ベクトルとしてのGF(2)疑似乱数生成器
     *
     * このクラスはGF(2)ベクトルとしての疑似乱数生成器を表す。また
     * GF(2)係数多項式をv個まとめたベクトルとしても扱うが、多項式特有の
     * 処理としては不定元倍(x をかける）ことのみをサポートする。これは
     * next_state()メソッドによって行われる。
     *
     * @tparam U 疑似乱数生成器の出力の型, 符号なし型でなければならない
     *\endjapanese
     *
     *\english
     * @brief GF(2) pseudo random number generator as a GF(2) vector.
     *
     * This class treats a GF(2) pseudo random number generator as a
     * GF(2) vector. And also the class treats the generator as a
     * vector of polynomials over GF(2).
     *
     * @tparam U type of output of pseudo random number
     * generator. Should be unsigned type.
     *\endenglish
     */
    template<typename U, typename SIMDGenerator>
    class simd_linear_generator_vector {
    public:

        /**
         *\japanese
         * コンストラクタ
         *
         * @param generator 均等分布次元を計算するGF(2)疑似乱数生成器
         *\endjapanese
         *
         * Constructor
         * @param generator GF(2) pseudo random number generator,
         * whose dimension of equi-distribution will be calculated.
         *\english
         *\endenglish
         */
        simd_linear_generator_vector<U, SIMDGenerator>(
            const SIMDGenerator& generator,
            SIMDInfo& info,
            bool lsb = false)
            {
            shared_ptr<SIMDGenerator> r(new SIMDGenerator(generator));
            rand = r;
            count = 0;
            zero = false;
            setZero(next);
            this->info = info;
            this->lsb = lsb;
#if defined(USE_SPECIAL)
            special = 0;
#endif
        }

        /**
         *\japanese
         * 標準基底コンストラクタ
         *
         * 標準基底を構成するベクトルのひとつを生成する。
         * 出力の特定のビットに1度だけ1を出力して、その後はずっと0を出力する
         * 疑似乱数生成器のコンストラクタ。例えば、1, 0, 0, 0, .... と出力する、
         * あるいは 8, 0, 0, 0, .... と出力するなど。
         *
         * @param generator GF(2)疑似乱数生成器
         * @param bit_pos 1 の位置
         *\endjapanese
         *
         *\english
         * Constructor for standard basis.
         *
         * Construct a vector which consists standard basis.
         * As a generator, this generates a number include only one bit once,
         * and then generate zero forever. For instance, it generates
         * 1, 0, 0, 0, ... or 8, 0, 0, 0, ...
         *\endenglish
         */
        simd_linear_generator_vector<U, SIMDGenerator>(
            const SIMDGenerator& generator,
            int bit_pos, SIMDInfo& info, bool lsb = false) {
            shared_ptr<SIMDGenerator> r(new SIMDGenerator(generator));
            rand = r;
            rand->setZero();
            count = 0;
            zero = false;
//            next = getOne<U>() << (bit_size<U>() - bit_pos - 1);
            setZero(next);
            setBitOfPos(&next, bit_size<U>() - bit_pos - 1, 1);
            this->info = info;
            this->lsb = lsb;
#if defined(USE_SPECIAL)
            special = 1;
#endif
#if defined(DEBUG)
            cout << "DEBUG:" << dec << bit_pos << ":"
                 << hex << next << endl;
#endif
        }

        void add(const simd_linear_generator_vector<U, SIMDGenerator>& src);
        void get_next(int bit_len);
        void next_state(int bit_len);
        void debug_print();

        /**
         *\japanese
         * GF(2)線形疑似乱数生成器
         *\endjapanese
         *
         *\english
         * A GF(2) linear pseudo random number generator
         *\endenglish
         */
        shared_ptr<SIMDGenerator> rand;
        /**
         *\japanese
         * next_state() が呼ばれた回数
         * これは多項式としてみた場合の次数に関係がある。
         *\endjapanese
         *
         *\english
         * number which shows counts of next_state() called.
         * As a polynomial, this is related to degree of polynomial.
         *\endenglish
         */
        int count;

        /**
         *\japanese
         * ゼロベクトルであるかどうかを示す。
         *\endjapanese
         *
         *\english
         * Shows if zero vector or not.
         *\endenglish
         */
        bool zero;

        /**
         *\japanese
         * 疑似乱数生成器の最新の出力（の上位vビット）または
         * 多項式の最高次の係数のなすベクトル
         *\endjapanese
         *
         *\english
         * latest output of pseudo random number generator, or, a
         * GF(2) vector consists of coefficients of highest degree
         * term of polynomial.
         *\endenglish
         */
        U next;
        SIMDInfo info;
        bool lsb;
#if defined(USE_SPECIAL)
        int special;
#endif
    };

    /**
     * @class AlgorithmEquidistribution
     *\japanese
     * @brief 疑似乱数生成器の均等分布次元を計算する
     *
     * PIS法(原瀬)によって疑似乱数生成器の出力の均等分布次元を計算する
     * アルゴリズム
     *
     * @tparam U 疑似乱数生成器の出力の型
     *\endjapanese
     *
     *\english
     * @brief Calculate dimension of equi-distribution of output of
     * pseudo random number generators.
     *
     * Algorithm that calculates dimension of equi-distribution of
     * output of pseudo random number generators using PIS
     * method[1](S. Harase).
     * @tparam type of output of pseudo random number generator.
     *\endenglish
     */
    template<typename U, typename SIMDGenerator>
    class AlgorithmSIMDEquidistribution {

        /**
         *\japanese
         * GF(2)ベクトルとしての疑似乱数生成器
         *\endjapanese
         *
         *\english
         * Pseudo random number generator as a vector.
         *\endenglish
         */
        typedef simd_linear_generator_vector<U, SIMDGenerator> linear_vec;
    public:

        /**
         *\japanese
         * コンストラクタ
         *
         * PIS法の特徴としてvビット精度均等分布次元を計算する際に、k(v+1)
         * の計算時の中間結果を利用してk(v)の計算の手間を省くことができる。この
         * クラスではその特徴を反映してk(v)のvをbit_len から1まで変化さ
         * せて一度に求めることができるようになっている。
         *
         * @param rand 均等分布次元計算可能な疑似乱数生成器
         * @param bit_length 均等分布次元を計算するMSBからのビット長, k(v)のv
         * の最初の値
         *\endjapanese
         *
         *\english
         * Constructor
         *
         * PIS method can calculates dimension of equi-distribution of
         * v-bit accuracy k(v) using intermediate result of
         * calculation of k(v+1). This class calculate multiple k(v)
         * varying \b v from \b bit_length to 1.
         *
         * @param rand pseudo random number generator
         * @param bit_length bit length from MSB to calculate dimension
         * of equi-distribution. This is first v of k(v).
         *\endenglish
         */
        AlgorithmSIMDEquidistribution(const SIMDGenerator& rand,
                                      int bit_length,
                                      SIMDInfo info,
                                      int maxbitsize,
                                      bool lsb = false) {
#if defined(DEBUG)
            cout << "AlgorithmSIMDEquidistribution constructer start" << endl;
            cout << "bit_length = " << dec << bit_length << endl;
#endif
            int bit_len = bit_length;
            int bit_size = bit_len * info.elementNo;
#if defined(DEBUG)
            cout << "bit_size = " << dec << bit_size << endl;
#endif
            this->info = info;
            size = bit_size + 1;
            basis = new linear_vec * [size];
            stateBitSize = maxbitsize;
            for (int i = 0; i < bit_size; i++) {
                basis[i] = new linear_vec(rand, i, info, lsb);
            }
            basis[bit_size] = new linear_vec(rand, info, lsb);
            basis[bit_size]->next_state(bit_len);
#if defined(DEBUG)
            cout << "zero = " << dec << basis[bit_size]->zero << endl;
            cout << "count = " << dec << basis[bit_size]->count << endl;
            cout << "AlgorithmSIMDEquidistribution constructer end" << endl;
#endif
        }

        /**
         *\japanese
         * デストラクタ
         *\endjapanese
         *\english
         * Destructor
         *\endenglish
         */
        ~AlgorithmSIMDEquidistribution() {
            for (int i = 0; i < size; i++) {
                delete basis[i];
            }
            delete[] basis;
        }

        int get_equidist(int bitLen);
    private:
        int get_equidist_main(int bit_len);

        /**
         *\japanese
         * 標準基底+1個のベクトルからなる配列。
         * basis という名前だが基底ではなく格子の生成集合という方が正しい。
         * k(v)を算出した時点では、ゼロベクトルとなっている
         * ベクトルを除けばv次元の格子の基底となっている。
         *\endjapanese
         *
         *\english
         * An array consists of standard basis plus one vector.
         * This array has a name \b basis, but this array is not
         * basis, this array should be called generating set of
         * lattice. But this array becomes basis of v dimensional
         * lattice when k(v) has calculated excluding zero vector.
         *\endenglish
         */
        linear_vec **basis;

        /**
         *\japanese
         * vビット精度均等分布次元の計算の v
         *\endjapanese
         *\english
         * \b v of dimension of equi-distribution with v-bit accuracy.
         *\endenglish
         */
        //int bit_len;

        /**
         *\japanese
         * 疑似乱数生成器の状態空間のビット数
         *\endjapanese
         *\english
         * Number of bits in internal state of pseudo random number
         * generator.
         *\endenglish
         */
        int stateBitSize;

        /**
         *\japanese
         * basis の配列の要素数
         *\endjapanese
         *
         *\english
         * Size of array \b basis.
         *\endenglish
         */
        SIMDInfo info;
        int size;
    };

#if defined(DEBUG)
    /**
     *\japanese
     * デバッグ出力
     *\endjapanese
     *\english
     * debug output
     *\endenglish
     */
    template<typename U, typename V>
    void simd_linear_generator_vector<U, V>::debug_print() {
        using namespace std;

        cout << "count = " << dec << count;
        cout << " zero = " << dec << zero;
        cout << " next = " << hex << next << endl;
    }
#else
    template<typename U, typename V>
    void simd_linear_generator_vector<U, V>::debug_print() {
    }
#endif

    /**
     *\japanese
     * vビット精度の均等分布次元を計算する。
     * v = \b bit_len から 1までの均等分布次元を計算して、\b veq[]
     * に入れる。返却値はv=1からbit_len までの均等分布次元の理論的上限との
     * 差の総和である。
     *
     * \warning AlgorithmEquidistribution をコンストラクトしてから、
     * get_all_equidist() または、get_equidist() のどちらか一方を１回し
     * か呼び出すことはできない。
     *
     * @param[out] veq v ビット精度の均等分布次元の配列
     * @return 実際のvビット精度の均等分布次元と理論的上限の差の総和
     *\endjapanese
     *
     *\english
     *
     * Calculate dimension of equi-distribution with v-bit accuracy.
     *
     * Calculate dimension of equi-distribution with v-bit accuracy
     * k(v) for v = \b bit_length to 1, and set them into an array \b
     * veq[].  The return value is sum of d(v)s, which are difference
     * between k(v) and theoretical upper bound at \b v.
     *
     * \warning Only one of get_all_equidist() or get_equidist() can
     * be called only one time. (Bad interface)
     *
     * @param[out] veq an array of k(v)
     * @return sum of d(v)s
     *
     *\endenglish
     */
    template<typename U, typename V>
    int AlgorithmSIMDEquidistribution<U, V>::get_equidist(int bitLen)
    {
        using namespace std;
        return get_equidist_main(bitLen);
    }

    /**
     *\japanese
     * ベクトルの加法
     * @param src このベクトルに足す相手のベクトル
     *\endjapanese
     *
     *\english
     * Vector addition
     * @param src source vector to be added to this vector
     *\endenglish
     */
    template<typename U, typename SIMDGenerator>
    void simd_linear_generator_vector<U, SIMDGenerator>::add(
        const simd_linear_generator_vector<U, SIMDGenerator>& src) {
        using namespace std;

        rand->add(*src.rand);
        next ^= src.next;
#if defined(USE_SPECIAL)
        special = 0;
#endif
        //count = src.count; // DEBUG OK?
    }

    /**
     * weight を考慮してnextを作る。prev もセットする。
     *
     */
    template<typename U, typename SIMDGenerator>
    void simd_linear_generator_vector<U, SIMDGenerator>::get_next(int bitSize) {
        using namespace std;
        U w = rand->generate();
#if defined(DEBUG) && 0
        cout << "w = " << w << endl;
#endif
        bitSize = bitSize * info.elementNo;
        setZero(next);
        int k = info.bitSize - 1;
        if (info.bitMode == 32) {
            if (lsb) {
                for (int i = 0; i < info.elementNo; i++) {
                    w.u[i] = reverse_bit(w.u[i]);
                }
            }
            for (int i = 0; i < info.elementNo; i++) {
                uint32_t mask = UINT32_C(0x80000000);
                for (int j = 0; j < bitSize; j += info.elementNo) {
                    if (w.u[i] & mask) {
                        setBitOfPos(&next, k, 1);
                    } else {
                        setBitOfPos(&next, k, 0);
                    }
                    k--;
                    mask = mask >> 1;
                }
            }
        } else {
            if (lsb) {
                for (int i = 0; i < info.elementNo; i++) {
                    w.u64[i] = reverse_bit(w.u64[i]);
                }
            }
            for (int i = 0; i < info.elementNo; i++) {
                uint64_t mask = UINT64_C(0x8000000000000000);
                for (int j = 0; j < bitSize; j += info.elementNo) {
                    if (w.u64[i] & mask) {
                        setBitOfPos(&next, k, 1);
                    } else {
                        setBitOfPos(&next, k, 0);
                    }
                    k--;
                    mask = mask >> 1;
                }
            }
        }
#if defined(DEBUG)
        if (!isZero(next)) {
            cout << "bitSize = " << dec << bitSize;
            cout << " info.elementNo = " << dec << info.elementNo << endl;
            cout << "get_next w = " << hex << w << endl;
            cout << "get_next next = " << next << endl;
        }
#endif
    }
    /**
     *\japanese
     * 疑似乱数生成器の状態遷移
     *
     * 多項式ベクトルとしてみると、すべての多項式を不定元倍する。
     *
     * @param bit_len MSB からの bit 長
     *\endjapanese
     *
     *\english
     * State transition of pseudo random number generator.
     *
     * As a vector of polynomial, multiply by an indeterminate.
     * @param bit_len bit length from MSB, or \b v of k(v) currently
     * calculating.
     *\endenglish
     */
    template<typename U, typename V>
    void simd_linear_generator_vector<U, V>::next_state(int bitLen) {
        using namespace std;

        if (zero) {
            return;
        }
        int zero_count = 0;
        get_next(bitLen);
        count++;
        while (isZero(next)) {
            zero_count++;
            if (zero_count > rand->bitSize() * 2) {
                zero = true;
                if (rand->isZero()) {
                    zero = true;
                }
                break;
            }
            get_next(bitLen);
            count++;
        }
    }

    /**
     *\japanese
     * PIS法によるvビット精度均等分布次元の計算のメインとなるメソッド
     *
     * @param v MSBからのビット長
     * @return v ビット精度均等分布次元
     *\endjapanese
     *
     *\english
     * Main method of calculation of dimension of equi-distribution
     * with v-bit accuracy.
     *
     * @param v v of k(v)
     * @return k(v)
     *\endenglish
     */
    template<typename U, typename SIMDGenerator>
    int AlgorithmSIMDEquidistribution<U, SIMDGenerator>::
    get_equidist_main(int v)
    {
        using namespace std;
        using namespace NTL;
#if defined(DEBUG)
        cout << "get_equidist_main start" << endl;
        cout << "v = " << dec << v << endl;
#endif
        int bitSize = v * info.elementNo;
        int pivot_index;
        int old_pivot = 0;

        pivot_index = calc_1pos(basis[bitSize]->next);
#if defined(DEBUG)
        cout << "get_equidist_main step 1" << endl;
#endif
        while (!basis[bitSize]->zero) {
#if defined(DEBUG)
            cout << "get_equidist_main step 2" << endl;
            cout << "basis[bitSize]->count = " << dec
                 << basis[bitSize]->count << endl;
#endif
#if defined(DEBUG) || 1
            //debug
#if 0
            if (pivot_index == 0) {
                cout << "pivot_index = " << dec << pivot_index << endl;
                cout << "next = " << hex << basis[bitSize]->next << endl;
            }
#endif
            if (pivot_index == -1) {
                cout << "pivot_index = " << dec << pivot_index << endl;
                cout << "zero = " << basis[bitSize]->zero << endl;
                cout << "next = " << hex << basis[bitSize]->next << endl;
                throw new std::logic_error("pivot error 0");
            }
            if (pivot_index >= bitSize) {
                cout << "pivot_index = " << dec << pivot_index << endl;
                cout << "bitSize = " << bitSize << endl;
                cout << "next = " << hex << basis[bitSize]->next << endl;
                throw new std::logic_error("pivot error 0.1");
            }
            if (pivot_index != calc_1pos(basis[pivot_index]->next)) {
                cerr << "pivot error 1" << endl;
                cerr << "pivot_index:" << dec << pivot_index << endl;
                cerr << "calc_1pos:" << dec
                     << calc_1pos(basis[pivot_index]->next) << endl;
                cerr << "next:" << hex << basis[pivot_index]->next << endl;
                for (int i = 0; i < bitSize; i++) {
                    cerr << dec << i << ":" << hex << basis[i]->next << endl;
                }
                throw new std::logic_error("pivot error 1");
            }
#endif
            // アルゴリズムとして、全部のcount を平均的に大きくしたい。
            // 従って count の小さい方を変化させたい
#if defined(USE_SPECIAL)
            // special は優先的に小さくしたい
            if (basis[bitSize]->count > basis[pivot_index]->count
                || (basis[bitSize]->count == basis[pivot_index]->count
                    && basis[pivot_index]->special != 0)) {
                swap(basis[bitSize], basis[pivot_index]);
            }
#else
            if (basis[bitSize]->count > basis[pivot_index]->count) {
                swap(basis[bitSize], basis[pivot_index]);
            }
#endif

#if defined(DEBUG)
            cout << "before add bitSize next = " << hex
                 << basis[bitSize]->next << endl;
            cout << "before add pivot   next = " << hex
                 << basis[pivot_index]->next << endl;
#endif
            basis[bitSize]->add(*basis[pivot_index]);
#if defined(DEBUG)
            cout << "after  add bitSize next = " << hex
                 << basis[bitSize]->next << endl;
#endif
            // add の結果 next の最後の1 は必ず 0 になる。
            // 全部0なら次の状態に進める。（内部でcount が大きくなる）
            if (isZero(basis[bitSize]->next)) {
                basis[bitSize]->next_state(v);
                pivot_index = calc_1pos(basis[bitSize]->next);
#if defined(DEBUG)
                cout << "zero" << endl;
                cout << "pivot_index = " << dec << pivot_index << endl;
                cout << "next = " << hex << basis[bitSize]->next << endl; // DBUG
                if (pivot_index >= bitSize) {
                    cout << "pivot_index = " << dec << pivot_index << endl;
                    cout << "bitSize = " << bitSize << endl;
                    cout << "next = " << hex << basis[bitSize]->next << endl;
                    throw new std::logic_error("pivot error 1.1");
                }
                if (basis[bitSize]->zero) {
                    cout << "loop exit condition" << endl;
                } else if (pivot_index == -1) {
                    cerr << "pivot error 1.1" << endl;
                    throw new std::logic_error("pivot error 1.2");
                }
#endif
            // 全部0でなければ、pivot_index は小さくなる。
            // pivot_index が 0 になれば最上位bit のみ1なので
            // 次の add で全部0になる。
            } else {
                old_pivot = pivot_index;
                pivot_index = calc_1pos(basis[bitSize]->next);
                if (pivot_index >= bitSize) {
                    cout << "pivot_index = " << dec << pivot_index << endl;
                    cout << "bitSize = " << bitSize << endl;
                    cout << "next = " << hex << basis[bitSize]->next << endl;
                    throw new std::logic_error("pivot error 2");
                }
                if (old_pivot <= pivot_index) {
                    cerr << "pivot error 2" << endl;
                    cerr << "old_pivot = " << dec << old_pivot << endl;
                    cerr << "pivot_index = " << dec << pivot_index << endl;
                    throw new std::logic_error("pivot error 2.1");
                }
            }
        }

        // 計算終了したので最長のベクトルを求める。（長いとはcountが少ないこと）
#if defined(DEBUG)
        for (int i = 0; i < bitSize; i++) {
            cout << dec << i << ": count = " << basis[i]->count << endl;
        }
#endif
        int min_count = basis[bitSize]->count;
//        int min_count = INT_MAX;
        for (int i = 0; i < bitSize; i++) {
            if (basis[i]->zero) {
                continue;
            }
#if defined(USE_SPECIAL)
            // trick
            if (basis[i]->count == 0 && basis[i]->special == 0) {
                basis[i]->count = 1;
            }
#endif
            if (min_count > basis[i]->count) {
                min_count = basis[i]->count;
            }
        }
        int result = min_count;
        if (result > stateBitSize / bitSize) {
//        if (result > stateBitSize / (bitSize / info.elementNo)) {
            cerr << "over theoretical bound" << endl;
            cout << "bitSize = " << dec << bitSize << endl;
            cout << "min_count = " << dec << min_count << endl;
            cout << "stateBitSize = " << dec << stateBitSize << endl;
            cout << "elementNo = " << dec << info.elementNo << endl;
            cout << "result = " << dec << result << endl;
            cerr << basis[0]->rand->getParamString() << endl;
#if 0
            for(int i = 0; i < size; i++) {
                basis[i]->debug_print();
            }
#endif
            throw new std::logic_error("over theoretical bound");
        }
        if (result <= 0) {
            cerr << "under bound" << endl;
            for (int i = 0; i < bitSize; i++) {
                cerr << dec << i << ": count = " << basis[i]->count << endl;
                if (basis[i]->count == 0) {
                    cerr << basis[i]->next << endl;
                }
            }
            throw new std::logic_error("under bound");
        }
#if defined(DEBUG)
        cout << "get_equidist_main end" << endl;
#endif
        return result;
    }

    // v を指定してそこだけ求める
    template<typename U, typename SIMDGenerator>
    int calc_SIMD_equidist(int v,
                           const SIMDGenerator& rand,
                           SIMDInfo& info,
                           int mexp,
                           bool lsb = false)
    {
        using namespace std;
        int state_inc;
        int weight_max = info.bitSize / 32;
        int state_max = weight_max;
        if (info.bitMode == 32) {
            state_inc = 1;
        } else {
            state_inc = 2;
        }
        int weight_dec = state_inc;
        int weight_start = weight_dec;
        if (info.fastMode) {
            weight_start = weight_max;
        }
        int veq = INT_MAX;
        int veq_weight = -1;
        // start mode 中で一番小さいもの
        for (int sm = 0; sm < state_max; sm += state_inc) {
            // weight mode 中で一番大きい物
            for (int wm = weight_start; wm <= weight_max; wm += weight_dec) {
#if 0
                cout << "start_mode = " << dec << sm;
                cout << " weight_mode = " << dec << wm;
#endif
                SIMDGenerator work = rand;
                work.setStartMode(sm);
                work.setWeightMode(wm);
                work.generate();

                //for (int v = 1; v <= bit_len; v++) {
                AlgorithmSIMDEquidistribution<U, SIMDGenerator>
                    ase(work, v, info, rand.bitSize(), lsb);
                int e = ase.get_equidist(v);
#if 0
                cout << " min_count = " << dec << e << endl;
#endif
                int e2 = e;
                if (info.bitMode == 32) {
                    e2 = e * info.elementNo - (info.elementNo - wm);
                } else { // 64
                    e2 = e * info.elementNo - (info.elementNo - wm / 2);
                }
                if (e2 > mexp / v) {
                    cerr << "over theoretical bound" << endl;
                    cout << "start_mode = " << dec << sm;
                    cout << " weight_mode = " << dec << wm;
                    cout << " mexp = " << dec << mexp;
                    cout << " e = " << dec << e2;
                    cout << " v = " << dec << v << endl;
                    throw new std::logic_error("over theoretical bound");
                }
#if 0
                cout << "\tk(" << dec << v << ") = " << dec
                     << e2 << endl;
#endif
                // max
                if (e2 > veq_weight) {
                    veq_weight = e2;
                }
                //}
            }
            if (veq > veq_weight) {
                veq = veq_weight;
            }
        }
        return veq;
    }

    template<typename U, typename SIMDGenerator>
    int calc_SIMD_equidistribution(const SIMDGenerator& rand,
                                   int veq[],
                                   int bit_len,
                                   SIMDInfo& info,
                                   int mexp,
                                   bool lsb = false)
    {
        int sum = 0;
        for (int i = 0; i < bit_len; i++) {
            veq[i] = calc_SIMD_equidist<U, SIMDGenerator>(i + 1, rand, info,
                                                          mexp, lsb);
            sum += mexp / (i + 1) - veq[i];
        }
        return sum;
    }

}
#endif // MTTOOLBOX_ALGORITHM_EQUIDISTRIBUTION_HPP
