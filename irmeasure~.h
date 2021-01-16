#ifndef __IRMEASURE__
#define __IRMEASURE__

#include "m_pd.h"
#include <math.h> /* serve? */
#include <stdio.h> /* serve? */
#include <stdlib.h> /* serve? */

#include "./lib/z_array.h"
#include "./lib/z_fft.h"
#include "./lib/z_core.h"
#include "./lib/z_utils.h"

#define PIRO_MAX_MEASURE_CHANS 128

typedef enum {
  NOISE_MODE_WHITE = 0,
  NOISE_MODE_BROWN = 1,
  NOISE_MODE_PINK = 2 
} t_noise_mode;

typedef struct _ess
{
  t_uint T;
  double K1, K2;

  double lo_f_act;
  double hi_f_act;

  double f1;
  double f2;

  double RT;
  double rf1, rf2;

  double fade_in;
  double fade_out;

  double sample_rate;
  double amp;

  t_uint num_amp_specifiers;

  double amp_specifier[34];
} t_ess;

typedef struct _mls
{
  t_uint32 feedback_mask;
  t_uint32 lfsr;

  t_uint32 T;
  t_uint32 order;

  double amp;
} t_mls;

typedef struct _xorshift
{
  t_uint32 w;
  t_uint32 x;
  t_uint32 y;
  t_uint32 z;
} t_xorshift;

typedef struct _noise_params
{
  double prev_output;
  double alpha;

  double alpha0;
  double alpha1;
  double alpha2;
  double alpha3;
  double alpha4;

  double b0;
  double b1;
  double b2;
  double b3;
  double b4;
  double b5;
  double b6;

  double amp;
  double sample_rate;

  double fade_in;
  double fade_out;
  double RT;

  t_noise_mode mode;
  t_xorshift gen;
  t_uint T;
} t_noise_params;

#endif /* __IRMEASURE__ */
