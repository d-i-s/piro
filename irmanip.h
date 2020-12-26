#ifndef __IRMANIP__
#define __IRMANIP__

#include "m_pd.h"
#include <math.h> /* */
#include <stdio.h> /* a cosa serve? */
#include <stdlib.h> /* a cosa serve? */

#include "./lib/z_array.h"
#include "./lib/z_fft.h"
#include "./lib/z_core.h"
#include "./lib/z_utils.h"
#include "./lib/z_matrix_math.h"

typedef enum {
  IRMANIP_ERROR = 0,
  IRMANIP_PHASE = 1,
  IRMANIP_INVERT = 2,
  IRMANIP_TRIM = 3,
  IRMANIP_AVERAGE = 4
} t_irmanip_operation;

#endif /* __IRMANIP__ */
