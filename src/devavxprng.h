#pragma once
#ifndef DEVAVXPRNG_H
#define  DEVAVXPRNG_H

#include "config.h"

#if HAVE_STDINT_H
#include <stdint.h>
#else
#pragma GCC error "do not have stdint.h"
#endif

#if HAVE_STDLIB_H
#include <stdlib.h>
#else
#pragma GCC error "do not have stdint.h"
#endif

#if HAVE_INTTYPES_H
#include <inttypes.h>
#else
#pragma GCC error "do not have inttypes.h"
#endif

#if HAVE_TIME_H
#include <time.h>
#else
#pragma GCC error "do not have time.h"
#endif

#if HAVE_UNISTD_H
#include <unistd.h>
#else
#pragma GCC error "do not have unistd.h"
#endif

#if HAVE_GETOPT_H
#include <getopt.h>
#else
#pragma GCC error "do not have getopt.h"
#endif

#if HAVE_LIMITS_H
#include <limits.h>
#else
#pragma GCC error "do not have limits.h"
#endif

#if HAVE_STRING_H
#include <string.h>
#else
#pragma GCC error "do not have string.h"
#endif

#include <iostream>
#include <iomanip>
#include <string>
#include <MTToolBox/util.hpp>


namespace MTToolBox {
    template<typename T>
        inline T make_msb_mask(int n)
    {
        return static_cast<T>(0);
    }
}
#endif //  DEVAVXPRNG_H
