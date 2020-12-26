#ifndef __U_MULTICHANCONVOLVE_STRUCT__
#define __U_MULTICHANCONVOLVE_STRUCT__

#include "m_types.h"
#include "m_convolve_latency_modes.h"
#include "m_convolve_errors.h"

#include "u_output_channel_convolve.h"


	
typedef struct _multi_channel_convolve
{
  t_output_channel_convolve *chan_convolvers[MAX_CHANS];
	
  t_uint num_in_chans;
  t_uint num_out_chans;
  long N2M;
	
  vFloat *in_temps[MAX_CHANS];
  vFloat *temp1;
  vFloat *temp2;
	
  t_memory temporary_memory;
	
} t_multi_channel_convolve;
	
#endif // __U_MULTICHANCONVOLVE_STRUCT__

		
/* N.B. pass 0 or less for out_chans to create a mono-to-mono parallel convolver */	
void multi_channel_convolve_free(t_multi_channel_convolve *x);
t_multi_channel_convolve *multi_channel_convolve_new(t_uint in_chans, t_uint out_chans, t_convolve_latency_mode latency_mode, t_uint max_length);

void multi_channel_convolve_clear(t_multi_channel_convolve *x, t_bool resize);
t_convolve_error multi_channel_convolve_resize(t_multi_channel_convolve *x, t_uint in_chan, t_uint out_chan, t_uint impulse_length);
t_convolve_error multi_channel_convolve_set(t_multi_channel_convolve *x, t_uint in_chan, t_uint out_chan, t_float *input, t_uint impulse_length, t_bool resize);
	
void multi_channel_convolve_process(t_multi_channel_convolve *x, t_sample **ins, t_sample **outs, t_sample *dry_gain, t_sample *wet_gain, t_uint vec_size, t_uint active_in_chans, t_uint active_out_chans);
