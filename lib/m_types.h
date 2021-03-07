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

#define FFTLOG2_TRIG_OFFSET ((t_uint) 3)
#define PASS_TRIG_OFFSET ((t_uint) 2)

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

#include <emmintrin.h>
typedef __m128  vFloat;
typedef __m128d vDouble;

#ifndef TARGET_INTEL
#if defined( __i386__ ) || defined ( __x86_64__ ) || defined(WIN_VERSION)
#define TARGET_INTEL
#endif
#endif

#ifdef TARGET_INTEL
#ifndef VECTOR_F64_128BIT
#define VECTOR_F64_128BIT
#endif
#endif

#ifdef TARGET_INTEL

#define float2vector                                    _mm_set1_ps

#define F32_VEC_MUL_OP					_mm_mul_ps
#define F32_VEC_ADD_OP					_mm_add_ps
#define F32_VEC_SUB_OP					_mm_sub_ps
#define F32_VEC_SHUFFLE					_mm_shuffle_ps

#define F32_VEC_ULOAD                                   _mm_loadu_ps
#define F32_VEC_USTORE                                  _mm_storeu_ps

#define F32_SHUFFLE_CONST(z, y, x, w)	((z<<6)|(y<<4)|(x<<2)|w)

#define double2vector                                    _mm_set1_pd

#define F64_VEC_MUL_OP  _mm_mul_pd
#define F64_VEC_ADD_OP  _mm_add_pd
#define F64_VEC_SUB_OP  _mm_sub_pd
#define F64_VEC_SHUFFLE _mm_shuffle_pd

#define F64_VEC_ULOAD   _mm_loadu_pd
#define F64_VEC_USTORE  _mm_storeu_pd

#define F64_SHUFFLE_CONST(y, x) ((y<<1)|x)

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

#if defined (NT)
#include <malloc.h>
#include <string.h> /* memset */

#define ALIGNED_MALLOC(x) _aligned_malloc(x, 16)
#define ALIGNED_REALLOC(x, nbytes) _aligned_realloc(x, nbytes, 16)
#define ALIGNED_FREE _aligned_free

#endif /* NT */

/* memory alignment */

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
#define CSET(a,b) cm_cset(a,b)
#define CEXP(a) cm_cexp(a)
#define CABS_SQ(a) cm_cabs_sq(a)
#define CONJ(a) cm_conj(a)
#define CADD(a, b) cm_cadd(a, b)
#define CSUB(a, b) cm_csub(a, b)
#define CMUL(a,b) cm_cmul(a,b)
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

  ret.real = a*c - b*d;
  ret.imag = a*d + b*c;

  return ret;
}

static t_complex_sample cm_cdiv(t_complex_sample in1, t_complex_sample in2)
{
  t_sample a = CREAL(in1);
  t_sample b = CIMAG(in1);
  t_sample c = CREAL(in2);
  t_sample d = CIMAG(in2);
  t_sample e = 1.0 / (c*c + d*d);
  
  t_sample real = (a*c + b*d) * e;
  t_sample imag = (b*c - a*d) * e;

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

  ret = cm_cset(a*cos(b),a*sin(b));

  return ret;
}

#endif /* __M_TYPES__ */
