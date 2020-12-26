#include "u_output_channel_convolve.h"


void output_channel_convolve_free(t_output_channel_convolve *x)
{
  t_uint i;
	
  if (!x)
    return;
	
  for (i = 0; i < x->num_in_chans; i++)
    zero_latency_convolve_free(x->convolvers[i]);
		 
  freebytes(x, sizeof(t_output_channel_convolve));
}


t_output_channel_convolve *output_channel_convolve_new(t_uint input_chans, t_uint max_length, t_convolve_latency_mode latency_mode)
{
  t_output_channel_convolve *x = (t_output_channel_convolve *)getbytes(sizeof(t_output_channel_convolve));
  t_uint i;
	
  if (!x)
    return 0;
	
  if (input_chans > MAX_CHANS)
    input_chans = MAX_CHANS;
		
  x->num_in_chans = 0;
	
  for (i = 0; i < input_chans; i++) {
    x->convolvers[i] = zero_latency_convolve_new(max_length, latency_mode);
		
    if (!x->convolvers[i]) {
      output_channel_convolve_free(x);
      return 0;
    }
		
    x->num_in_chans++;
  }
		
  return (x);
}


t_convolve_error output_channel_convolve_resize(t_output_channel_convolve *x, t_uint in_chan, t_uint impulse_length)
{
  if (in_chan < x->num_in_chans) {
    if (!zero_latency_convolve_resize(x->convolvers[in_chan], impulse_length, false))
      return CONVOLVE_ERR_MEM_UNAVAILABLE;
  }
  else
    return CONVOLVE_ERR_IN_CHAN_OUT_OF_RANGE;
  
  return CONVOLVE_ERR_NONE;
}


t_convolve_error output_channel_convolve_set(t_output_channel_convolve *x, t_uint in_chan, t_float *input, t_uint impulse_length, t_bool resize)
{
  if (in_chan < x->num_in_chans)
    return zero_latency_convolve_set(x->convolvers[in_chan], input, impulse_length, resize);
  else
    return CONVOLVE_ERR_IN_CHAN_OUT_OF_RANGE;
}


void output_channel_convolve_process(t_output_channel_convolve *x, vFloat **ins, t_sample *out, vFloat *temp1, vFloat *temp2, t_uint vec_size, t_uint active_in_chans)
{
  t_float32 *out_temp = (t_float32 *)temp2;
	
  t_uint i, j;
	
  /* Zero out temp */
  for (j = 0; j < vec_size; j++)
    out_temp[j] = 0.f;

  /* Convolve */
  for (i = 0; i < x->num_in_chans && i < active_in_chans; i++)
    zero_latency_convolve_process(x->convolvers[i], ins[i], temp1, temp2, vec_size);
	
  /* Copy output */
  for (j = 0; j < vec_size; j++)
    out[j] = (t_sample)out_temp[j];
}
