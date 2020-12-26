#ifndef __U_TIMECONVOLVE_STRUCT__
#define __U_TIMECONVOLVE_STRUCT__

#include "m_convolve_errors.h"
#include "m_types.h"

#include <string.h>

typedef struct _time_domain_convolve
{
  /* Internal buffers */
	
  t_float32 *impulse_buffer;
  t_float32 *input_buffer;
	
  t_uint input_position; 
  t_uint impulse_length;
	
  t_uint offset; 
  t_uint length;
	
} t_time_domain_convolve;

#endif // __U_TIMECONVOLVE_STRUCT__


void time_domain_convolve_free(t_time_domain_convolve *x);
t_time_domain_convolve *time_domain_convolve_new(t_uint offset, t_uint length);

t_convolve_error time_domain_convolve_length_set(t_time_domain_convolve *x, t_uint length);
void time_domain_convolve_offset_set(t_time_domain_convolve *x, t_uint offset);

t_convolve_error time_domain_convolve_set(t_time_domain_convolve *x, t_float *input, t_uint impulse_length);

void time_domain_convolve_process_scalar(t_time_domain_convolve *x, t_float32 *in, t_float32 *out, t_uint vec_size);
void time_domain_convolve_process(t_time_domain_convolve *x, t_float32 *in, t_float32 *out, t_uint vec_size);
