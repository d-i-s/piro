#include "m_types.h"

#include "u_zero_latency_convolve.h"
#include "m_convolve_errors.h"

#ifndef __OUTPUTCHANCONVOLVE_STRUCT__
#define __OUTPUTCHANCONVOLVE_STRUCT__

#define MAX_CHANS 96

typedef struct _output_channel_convolve
{
  t_zero_latency_convolve *convolvers[MAX_CHANS];
	
  t_uint num_in_chans;
		
} t_output_channel_convolve;

#endif //__OUTPUTCHANCONVOLVE_STRUCT__

void output_channel_convolve_free(t_output_channel_convolve *x);
t_output_channel_convolve *output_channel_convolve_new(t_uint input_chans, t_uint max_length, t_convolve_latency_mode latency_mode);

t_convolve_error output_channel_convolve_resize(t_output_channel_convolve *x, t_uint in_chan, t_uint impulse_length);
t_convolve_error output_channel_convolve_set(t_output_channel_convolve *x, t_uint in_chan, t_float *input, t_uint impulse_length, t_bool resize);

void output_channel_convolve_process(t_output_channel_convolve *x, vFloat **ins, t_sample *out, vFloat *temp1, vFloat *temp2, t_uint vec_size, t_uint active_in_chans);
