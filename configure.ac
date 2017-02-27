#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([DevAVXPRNG], [0.1.0], [saito@manieth.com])
AC_CONFIG_SRCDIR([include/config.h.in])
AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL

AX_GCC_ARCHFLAG

# Open MPI compiler
AC_ARG_ENABLE([mpi],
        AS_HELP_STRING([--enable-mpi], [make open mpi search executables]))
AM_CONDITIONAL([MPI], [test x$enable_mpi = xyes])

if test !AS_VAR_TEST_SET([MPICXX]); then
  AS_VAR_SET([MPICXX], [`which mpicxx`])
fi
AC_SUBST(MPICXX)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([inttypes.h stdint.h stdlib.h string.h memory.h \
x86intrin.h immintrin.h zmmintrin.h cpuid.h time.h getopt.h limits.h \
string.h])
#AX_EXT

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT64_T
AC_TYPE_UINT64_T

# keep this order
AC_CHECK_LIB(gf2x, gf2x_mul)
AC_CHECK_LIB(gmp, __gmpz_init)
AC_CHECK_LIB(ntl, main,[],[AC_MSG_ERROR(Need ntl)])
AC_CHECK_LIB(MTToolBox, get_mttoolbox_version,[],[AC_MSG_ERROR(Need MTToolBox)])

AX_CXX_COMPILE_STDCXX_0X

AX_CXX_COMPILE_STDCXX_14(noext, optional) # keep this order
AX_CXX_COMPILE_STDCXX_11(noext, optional)

# Checks for library functions.
#AC_CHECK_FUNCS([strtoull _may_i_use_cpu_feature __builtin_cpu_supports])
AC_MSG_CHECKING(for _may_i_use_cpu_feature)
AC_DEFINE(HAVE_MAY_I_USE_CPU_FEATURE, [0],
[Define to 1 if you have _may_i_use_cpu_feature function])
AC_CACHE_VAL(ac_cv_may_i_use_cpu_feature,
AC_TRY_COMPILE([],
 [_may_i_use_cpu_feature(0)],
 ac_cv_may_i_use_cpu_feature=yes,
 ac_cv_may_i_use_cpu_feature=no))
AC_MSG_RESULT($ac_cv_may_i_use_cpu_feature)
if test $ac_cv_may_i_use_cpu_feature = yes; then
  AC_DEFINE(HAVE_MAY_I_USE_CPU_FEATURE, [1],
  [Define to 1 if you have _may_i_use_cpu_feature function])
fi

AC_MSG_CHECKING(for __builtin_cpu_supports)
AC_DEFINE(HAVE_BUILTIN_CPU_SUPPORTS, [0],
[Define to 1 if you have __builtin_cpu_supports function])
AC_CACHE_VAL(ac_cv_builtin_cpu_supports,
AC_TRY_COMPILE([],
 [__builtin_cpu_supports("mmx")],
 ac_cv_builtin_cpu_supports=yes,
 ac_cv_builtin_cpu_supports=no))
AC_MSG_RESULT($ac_cv_builtin_cpu_supports)
if test $ac_cv_builtin_cpu_supports = yes; then
  AC_DEFINE(HAVE_BUILTIN_CPU_SUPPORTS, [1],
  [Define to 1 if you have __builtin_cpu_supports function])
fi

AC_MSG_CHECKING(for __builtin_parityll)
AC_DEFINE(HAVE_BUILTIN_PARITYLL, [0],
[Define to 1 if you have __builtin_parityll function])
AC_CACHE_VAL(ac_cv_builtin_parityll,
AC_TRY_COMPILE([],
 [__builtin_parityll(0)],
 ac_cv_builtin_parityll=yes,
 ac_cv_builtin_parityll=no))
AC_MSG_RESULT($ac_cv_builtin_parityll)
if test $ac_cv_builtin_parityll = yes; then
  AC_DEFINE(HAVE_BUILTIN_PARITYLL, [1],
  [Define to 1 if you have __builtin_parityll function])
fi

AX_CHECK_COMPILE_FLAG([-mavx512f], [AX_APPEND_FLAG([-mavx512f], [CXXFLAGS])],
 [], [],
 [AC_LANG_PROGRAM([#include "zmmintrin.h"], [__m512i a, b;_mm512_xor_si512(a,b)])])
AX_CHECK_COMPILE_FLAG([-mavx2], [AX_APPEND_FLAG([-mavx2], [CXXFLAGS])],
 [], [],
 [AC_LANG_PROGRAM([#include "immintrin.h"], [__m256i a,b;_mm256_zeroupper();_mm256_xor_si256(a,b)])])
AX_CHECK_COMPILE_FLAG([-msse4.2], [AX_APPEND_FLAG([-msse4.2], [CXXFLAGS])],
 [], [],
 [AC_LANG_PROGRAM([#include "nmmintrin.h"], [_mm_crc32_u16(0,0)])])
AX_CHECK_COMPILE_FLAG([-msse4.1], [AX_APPEND_FLAG([-msse4.1], [CXXFLAGS])],
 [], [],
 [AC_LANG_PROGRAM([#include "smmintrin.h"], [__m128i a,b;_mm_mul_epi32(a,b)])])
AX_CHECK_COMPILE_FLAG([-mssse3], [AX_APPEND_FLAG([-mssse3], [CXXFLAGS])],
 [], [],
 [AC_LANG_PROGRAM([#include "tmmintrin.h"], [__m128i a;_mm_abs_epi16(a)])])
AX_CHECK_COMPILE_FLAG([-msse2], [AX_APPEND_FLAG([-msse2], [CXXFLAGS])],
 [], [],
 [AC_LANG_PROGRAM([#include "emmintrin.h"], [__m128i a,b;_mm_and_si128(a,b)])])

AC_DEFINE(HAVE_AVX512F, [0], [Define to 1 if compiler can use -mavx512f])
AC_DEFINE(HAVE_AVX2, [0], [Define to 1 if compiler can use -mavx2])
AC_DEFINE(HAVE_SSE4_2, [0], [Define to 1 if compiler can use -msse4.2])
AC_DEFINE(HAVE_SSE4_1, [0], [Define to 1 if compiler can use -msse4.1])
AC_DEFINE(HAVE_SSSE3, [0], [Define to 1 if compiler can use -mssse3])
AC_DEFINE(HAVE_SSE2, [0], [Define to 1 if compiler can use -msse2])
for flag in $CXXFLAGS; do
  if test $flag = "-mavx512f"; then
    AC_DEFINE(HAVE_AVX512F, [1], [Define to 1 if compiler can use -mavx512f])
  fi
  if test $flag = "-mavx2"; then
    AC_DEFINE(HAVE_AVX2, [1], [Define to 1 if compiler can use -mavx2])
  fi
  if test $flag = "-msse4.2"; then
    AC_DEFINE(HAVE_SSE4_2, [1], [Define to 1 if compiler can use -msse4.2])
  fi
  if test $flag = "-msse4.2"; then
    AC_DEFINE(HAVE_SSE4_1, [1], [Define to 1 if compiler can use -msse4.1])
  fi
  if test $flag = "-mssse3"; then
    AC_DEFINE(HAVE_SSSE3, [1], [Define to 1 if compiler can use -mssse3])
  fi
  if test $flag = "-msse2"; then
    AC_DEFINE(HAVE_SSE2, [1], [Define to 1 if compiler can use -msse2])
  fi
done

AC_DEFINE(HAVE_STD_SP, [0], [Define to 1 if you have std::shared_ptr])
AC_DEFINE(HAVE_STD_TR1_SP, [0], [Define to 1 if you have std::tr1::shared_ptr])
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_MSG_CHECKING(for std::shared_ptr)
AC_TRY_COMPILE([#include <memory>],
[std::shared_ptr<void> sp;return 0;],
ac_cv_cxx_have_std_sp=yes, ac_cv_cxx_have_std_sp=no)
AC_MSG_RESULT($ac_cv_cxx_have_std_sp)
AC_MSG_CHECKING(for std::tr1::shared_ptr)
AC_TRY_COMPILE([#include <tr1/memory>],
[std::tr1::shared_ptr<void> sp;return 0;],
ac_cv_cxx_have_std_tr1_sp=yes, ac_cv_cxx_have_std_tr1_sp=no)
AC_MSG_RESULT($ac_cv_cxx_have_std_tr1_sp)
AC_LANG_RESTORE
if test $ac_cv_cxx_have_std_sp = yes; then
  AC_DEFINE(HAVE_STD_SP, [1],
  [Define to 1 if you have std::shared_ptr])
fi
if test $ac_cv_cxx_have_std_tr1_sp = yes; then
  AC_DEFINE(HAVE_STD_TR1_SP, [1],
  [Define to 1 if you have std::tr1::shared_ptr])
fi


AC_CONFIG_FILES([Makefile src/Makefile mpisrc/Makefile])
AC_OUTPUT
