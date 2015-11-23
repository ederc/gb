#ifndef GB_GB_CONFIG_H
#define GB_GB_CONFIG_H

#include <stdint.h>
#include <stdio.h>
#include <limits.h>

/* Defines determined during configuration of gb. */
#define __GB_HAVE_MM_MALLOC		1
#define __GB_HAVE_POSIX_MEMALIGN	1
#define __GB_HAVE_SSE2		1
#define __GB_HAVE_OPENMP		1
#define __GB_CPU_L1_CACHE		
#define __GB_CPU_L2_CACHE		
#define __GB_CPU_L3_CACHE		
#define __GB_DEBUG_DUMP		(0 || @GB_DEBUG_MZD@)
#define __GB_DEBUG_MZD		@GB_DEBUG_MZD@
#define __GB_HAVE_LIBPNG              1

#define __GB_CC                       "gcc"
#define __GB_CFLAGS                   " -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -maes -mavx -mfma -mavx2 -fopenmp -g -O2 -DNDEBUG"
#define __GB_SIMD_CFLAGS              " -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -maes -mavx -mfma -mavx2"
#define __GB_OPENMP_CFLAGS            "-fopenmp"

/* Helper macros. */
#define __GB_USE_MM_MALLOC		(__GB_HAVE_MM_MALLOC && __GB_HAVE_SSE2)
#define __GB_USE_POSIX_MEMALIGN	(__GB_HAVE_POSIX_MEMALIGN && __GB_HAVE_SSE2)
#define __GB_DD_QUIET			(@GB_DEBUG_MZD@ && !0)

#define __GB_LOOP_UNROLL_SMALL  16
#define __GB_LOOP_UNROLL_BIG    64

#define __GB_ROUND_DOWN(x, s) ((x) & ~((s)-1))

// Note: It must hold that __GB_SIMD_BLOCK_SIZE % __GB_SIMD_INNER_SIZE == 0
#define __GB_COLUMN_B  1
#define __GB_SIMD_INNER_SIZE  4
#define __GB_SIMD_BLOCK_SIZE  256
#define __GB_SIMD_INNER_BLOCKS_PER_ROW  ((__GB_SIMD_BLOCK_SIZE) / (__GB_SIMD_INNER_SIZE))
#define __GB_SIMD_BLOCK_SIZE_RECT ((__GB_SIMD_BLOCK_SIZE) * (__GB_SIMD_BLOCK_SIZE))
#define __GB_SIMD_BLOCK_SIZE_DIAG (((__GB_SIMD_BLOCK_SIZE) *((__GB_SIMD_BLOCK_SIZE)-1)/2) + \
  (__GB_SIMD_BLOCK_SIZE))
#define __GB_NROWS_MULTILINE  2
#define __GB_HYBRID_THRESHOLD 0.5f
#define __GB_DENSITY_THRESHOLD 10.0f

#define __GB_WORDSIZE __WORDSIZE

static const uint32_t __GB_MINUS_ONE_32 = (uint32_t) -1;
static const uint16_t __GB_MINUS_ONE_16 = (uint16_t) -1;
static const uint8_t  __GB_MINUS_ONE_8  = (uint8_t)  -1;
#endif
