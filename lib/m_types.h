/*
 ARM support added by charlesneimog 26/04/13
*/

#ifndef __M_TYPES__
#define __M_TYPES__

#include "m_pd.h"
#include <math.h> /* sin, cos, M_PI, ... */

#define DEBUG32BIT (1) /* questo definisce il debug per i 32 bit */
#define DEBUG32BIT2 (1)

/* PD_FLOATSIZE = 32 definisce t_sample float o double */
/* PD_LONGINTTYPE long o long long --> t_int */
/* PD_FLOATTYPE float o double --> t_float, t_floatarg, t_sample */
/* PD_FLOATUINTTYPE unsigned int o unsigned long */

#ifdef PD_FLOAT_PRECISION
#define PD_FLOATSIZE 64
#endif /* PD_FLOAT_PRECISION */

#if PD_FLOATSIZE == 64
typedef unsigned long long t_uint;
#elif PD_FLOATSIZE == 32
typedef unsigned long t_uint;
#define PD_SAMPLE_PRECISION
#endif /* PD_FLOATSIZE */

#define SQRT_2_2 0.70710678118654752440084436210484904

#define UNSIGNED_INT32_TO_NORM_DOUBLE 2.32830643653869628906e-10

typedef unsigned int t_uint32;

#define FFTLOG2_TRIG_OFFSET ((t_uint)3)
#define PASS_TRIG_OFFSET ((t_uint)2)

typedef int t_bool;
#undef true
#define true 1
#undef false
#define false 0

typedef float t_float32; /* for partition convolver */

/* FFT data structure */

typedef struct _FFT_Split
{
  t_sample *realp;
  t_sample *imagp;
} FFT_Split;

typedef struct _FFT_Setup
{
  t_uint max_fft_log2;
  FFT_Split tables[28];
} FFT_Setup;

typedef struct _FFT_Split32
{
  t_float32 *realp;
  t_float32 *imagp;
} FFT_Split32;

typedef struct _FFT_Setup32
{
  t_int max_fft_log2; // t_uint
  FFT_Split32 tables[28];
} FFT_Setup32;

/* SIMD support */

#define ENABLE_SIMD_SUPPORT (0)

/* ------------------------------------------------------------------
   ADDITIONS FOR APPLE SILICON / ARM (non-destructive to existing code)
   ------------------------------------------------------------------ */

/* Detect Apple ARM64 (macOS on Apple Silicon) */
#if defined(__APPLE__) && (defined(__arm64__) || defined(__aarch64__))
#define TARGET_APPLE_ARM64 1
#endif

/* Generic ARM NEON detection */
#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(TARGET_APPLE_ARM64)
#define TARGET_ARM_NEON 1
#endif

/* Keep original TARGET_INTEL detection (unchanged) */
#ifndef TARGET_INTEL
#if defined(__i386__) || defined(__x86_64__) || defined(WIN_VERSION)
#define TARGET_INTEL
#endif
#endif

/* ------------- ORIGINAL INTEL SSE2 SECTION (kept, only wrapped) -------------
   (We only wrap the unconditional include so non-Intel targets still compile) */
#ifdef TARGET_INTEL
#include <emmintrin.h>
typedef __m128 vFloat;
typedef __m128d vDouble;

#ifndef VECTOR_F64_128BIT
#define VECTOR_F64_128BIT
#endif
#endif /* TARGET_INTEL */

/* If NOT Intel but ARM NEON is available, define vector types and macros */
#if !defined(TARGET_INTEL) && defined(TARGET_ARM_NEON)
#include <arm_neon.h>

/* Map SSE-like abstractions to NEON for float */
typedef float32x4_t vFloat;
/* Double precision NEON (only guaranteed on AArch64) */
#if defined(__aarch64__)
typedef float64x2_t vDouble;
#define HAVE_NEON_F64 1
#else
/* Fallback struct if 64-bit NEON not present (32-bit ARM) */
typedef struct
{
  double d[2];
} vDouble;
#endif

/* Float broadcast */
#define float2vector vdupq_n_f32
/* Basic float ops */
#define F32_VEC_MUL_OP vmulq_f32
#define F32_VEC_ADD_OP vaddq_f32
#define F32_VEC_SUB_OP vsubq_f32
#define F32_VEC_ULOAD(p) vld1q_f32((const float *)(p))
#define F32_VEC_USTORE(p, v) vst1q_f32((float *)(p), (v))

/* Limited shuffle emulation:
   SSE _mm_shuffle_ps(a,b,imm8) selects lanes from (a,b).
   We approximate using vext (concatenate + extract). This
   handles many common patterns (rotations / simple mixes). */
#define F32_VEC_SHUFFLE(a, b, imm8) vextq_f32((a), (b), ((imm8) & 0x3))

#define F32_SHUFFLE_CONST(z, y, x, w) ((z << 6) | (y << 4) | (x << 2) | w) /* kept for compatibility */

/* Double broadcast / ops (only if real NEON F64) */
#if defined(HAVE_NEON_F64)
#define double2vector vdupq_n_f64
#define F64_VEC_MUL_OP vmulq_f64
#define F64_VEC_ADD_OP vaddq_f64
#define F64_VEC_SUB_OP vsubq_f64
#define F64_VEC_ULOAD(p) vld1q_f64((const double *)(p))
#define F64_VEC_USTORE(p, v) vst1q_f64((double *)(p), (v))
#define F64_VEC_SHUFFLE(a, b, imm) vextq_f64((a), (b), ((imm) & 0x1))
#else
/* Scalar fallbacks if double NEON not available */
static inline vDouble double2vector(double x)
{
  vDouble r;
  r.d[0] = x;
  r.d[1] = x;
  return r;
}
static inline vDouble F64_VEC_MUL_OP(vDouble a, vDouble b)
{
  vDouble r;
  r.d[0] = a.d[0] * b.d[0];
  r.d[1] = a.d[1] * b.d[1];
  return r;
}
static inline vDouble F64_VEC_ADD_OP(vDouble a, vDouble b)
{
  vDouble r;
  r.d[0] = a.d[0] + b.d[0];
  r.d[1] = a.d[1] + b.d[1];
  return r;
}
static inline vDouble F64_VEC_SUB_OP(vDouble a, vDouble b)
{
  vDouble r;
  r.d[0] = a.d[0] - b.d[0];
  r.d[1] = a.d[1] - b.d[1];
  return r;
}
static inline vDouble F64_VEC_ULOAD(const double *p)
{
  vDouble r;
  r.d[0] = p[0];
  r.d[1] = p[1];
  return r;
}
static inline void F64_VEC_USTORE(double *p, vDouble v)
{
  p[0] = v.d[0];
  p[1] = v.d[1];
}
#define F64_VEC_SHUFFLE(a, b, imm) (a) /* no-op placeholder */
#endif                                 /* HAVE_NEON_F64 */

#define F64_SHUFFLE_CONST(y, x) ((y << 1) | x)

#endif /* !TARGET_INTEL && TARGET_ARM_NEON */

/* If both architectures somehow defined, prefer Intel originals (no changes). */

/* Keep original Intel macro block (unchanged) but guard it so we don't
   redefine when we already created ARM versions. */
#ifdef TARGET_INTEL

#define float2vector _mm_set1_ps

#define F32_VEC_MUL_OP _mm_mul_ps
#define F32_VEC_ADD_OP _mm_add_ps
#define F32_VEC_SUB_OP _mm_sub_ps
#define F32_VEC_SHUFFLE _mm_shuffle_ps

#define F32_VEC_ULOAD _mm_loadu_ps
#define F32_VEC_USTORE _mm_storeu_ps

#define F32_SHUFFLE_CONST(z, y, x, w) ((z << 6) | (y << 4) | (x << 2) | w)

#define double2vector _mm_set1_pd

#define F64_VEC_MUL_OP _mm_mul_pd
#define F64_VEC_ADD_OP _mm_add_pd
#define F64_VEC_SUB_OP _mm_sub_pd
#define F64_VEC_SHUFFLE _mm_shuffle_pd

#define F64_VEC_ULOAD _mm_loadu_pd
#define F64_VEC_USTORE _mm_storeu_pd

#define F64_SHUFFLE_CONST(y, x) ((y << 1) | x)

#endif /* TARGET_INTEL */

/*

__arm__ o __ARM_ARCH_7__

__m128 float32x4_t
_mm_loadu_ps vld1q_f32
_mm_storeu_ps vst1q_f32
_mm_add_ps vaddq_f32
_mm_sub_ps vsubq_f32
_mm_mul_ps vmulq_f32
_mm_set1_ps vdupq_n_f32


__128d float64x2_t
vaddq_f64
vsubq_f64
vmulq_f64
vld1q_f64
vst1q_f64
vdupq_n_f64


 */

/* static int SSE2_check(); */

#if defined(NT)
#include <malloc.h>
#include <string.h> /* memset */

#define ALIGNED_MALLOC(x) _aligned_malloc(x, 16)
#define ALIGNED_REALLOC(x, nbytes) _aligned_realloc(x, nbytes, 16)
#define ALIGNED_FREE _aligned_free

#else
/* POSIX (Linux, macOS including Apple ARM) aligned allocation helpers (added) */
#include <stdlib.h>
#include <string.h>
static inline void *aligned_malloc(size_t size, size_t alignment)
{
  void *ptr = NULL;
#if defined(_ISOC11_SOURCE)
  /* Ensure size is a multiple of alignment for aligned_alloc */
  size_t padded = ((size + alignment - 1) / alignment) * alignment;
  ptr = aligned_alloc(alignment, padded);
#elif defined(_MSC_VER)
  ptr = _aligned_malloc(size, alignment);
#else
  if (posix_memalign(&ptr, alignment, size) != 0)
    ptr = NULL;
#endif
  return ptr;
}
#define ALIGNED_MALLOC(x) aligned_malloc((x), 16)
#define ALIGNED_REALLOC(x, nbytes) realloc((x), (nbytes)) /* alignment may be lost */
#define ALIGNED_FREE free
#endif /* NT */

/* memory alignment (original declarations kept) */

void *aligned_getbytes(size_t nbytes);
void *aligned_resizebytes(void *old, size_t oldsize, size_t newsize);
void aligned_freebytes(void *fatso, size_t nbytes);

/* complex math */

typedef struct _complex_sample
{
  t_sample real;
  t_sample imag;
} t_complex_sample;

#define CREAL(a) (a).real
#define CIMAG(a) (a).imag
#define CSET(a, b) cm_cset(a, b)
#define CEXP(a) cm_cexp(a)
#define CABS_SQ(a) cm_cabs_sq(a)
#define CONJ(a) cm_conj(a)
#define CADD(a, b) cm_cadd(a, b)
#define CSUB(a, b) cm_csub(a, b)
#define CMUL(a, b) cm_cmul(a, b)
#define CDIV(a, b) cm_cdiv(a, b)

static t_complex_sample cm_cset(t_sample a, t_sample b)
{
  t_complex_sample ret;

  ret.real = a;
  ret.imag = b;

  return ret;
}

static t_complex_sample cm_cadd(t_complex_sample in1, t_complex_sample in2)
{
  t_sample a = CREAL(in1);
  t_sample b = CIMAG(in1);
  t_sample c = CREAL(in2);
  t_sample d = CIMAG(in2);

  return CSET(a + c, b + d);
}

static t_complex_sample cm_csub(t_complex_sample in1, t_complex_sample in2)
{
  t_sample a = CREAL(in1);
  t_sample b = CIMAG(in1);
  t_sample c = CREAL(in2);
  t_sample d = CIMAG(in2);

  return CSET(a - c, b - d);
}

static t_complex_sample cm_cmul(t_complex_sample in1, t_complex_sample in2)
{
  t_complex_sample ret;

  t_sample a = in1.real;
  t_sample b = in1.imag;
  t_sample c = in2.real;
  t_sample d = in2.imag;

  ret.real = a * c - b * d;
  ret.imag = a * d + b * c;

  return ret;
}

static t_complex_sample cm_cdiv(t_complex_sample in1, t_complex_sample in2)
{
  t_sample a = CREAL(in1);
  t_sample b = CIMAG(in1);
  t_sample c = CREAL(in2);
  t_sample d = CIMAG(in2);
  t_sample e = 1.0 / (c * c + d * d);

  t_sample real = (a * c + b * d) * e;
  t_sample imag = (b * c - a * d) * e;

  return CSET(real, imag);
}

static t_sample cm_cabs_sq(t_complex_sample in)
{
  t_sample a = CREAL(in);
  t_sample b = CIMAG(in);

  return (a * a) + (b * b);
}

static t_complex_sample cm_conj(t_complex_sample in)
{
  t_complex_sample ret;

  ret.real = in.real;
  ret.imag = -in.imag;

  return ret;
}

static t_complex_sample cm_cexp(t_complex_sample in)
{
  t_complex_sample ret;

  t_sample a = exp(in.real);
  t_sample b = in.imag;

  ret = cm_cset(a * cos(b), a * sin(b));

  return ret;
}

#endif /* __M_TYPES__ */