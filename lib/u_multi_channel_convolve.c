#include "u_multi_channel_convolve.h"


void multi_channel_convolve_free(t_multi_channel_convolve *x)
{
  t_uint i;
	
  if (!x)
    return;
	
  for (i = 0; i < x->num_out_chans; i++)
    output_channel_convolve_free(x->chan_convolvers[i]);
	
  /* aligned_freebytes(x->in_temps[0], 1); */

  /* free temporary memory */
  clear_memory(&x->temporary_memory);

  freebytes(x, sizeof(t_multi_channel_convolve));
}


void multi_channel_convolve_temp_setup(t_multi_channel_convolve *x, void *mem_pointer, t_uint max_frame_size)
{
  t_uint num_in_chans = x->num_in_chans;
  t_uint i;
	
  x->in_temps[0] = mem_pointer;
	
  for (i = 1; i < num_in_chans; i++)
    x->in_temps[i] = x->in_temps[0] + (i * (max_frame_size >> 2));
	
  x->temp1 = x->in_temps[num_in_chans - 1] + (max_frame_size >> 2);
  x->temp2 = x->temp1 + (max_frame_size >> 2);
}


t_multi_channel_convolve *multi_channel_convolve_new(t_uint num_in_chans, t_uint num_out_chans, t_convolve_latency_mode latency_mode, t_uint max_length)
{
  t_multi_channel_convolve *x = (t_multi_channel_convolve *)getbytes(sizeof(t_multi_channel_convolve));
  long N2M;
  t_uint i;
	
  if (!x)
    return 0;

  num_in_chans = num_in_chans < 1 ? 1 : num_in_chans;
  num_in_chans = num_in_chans > MAX_CHANS ? MAX_CHANS : num_in_chans;
  num_out_chans = num_out_chans > MAX_CHANS ? MAX_CHANS : num_out_chans;

  x->N2M = N2M = num_out_chans > 0;
  x->num_in_chans = 0;
  x->num_out_chans = 0;
  x->in_temps[0] = NULL;
	
  num_out_chans = !N2M ? num_in_chans : num_out_chans;
	
  for (i = 0; i < num_out_chans; i++) {
    x->chan_convolvers[i] = output_channel_convolve_new(N2M ? num_in_chans : 1, max_length, latency_mode);
		
    if (!x->chan_convolvers[i]) {
      multi_channel_convolve_free(x);
      return 0;
    }
		
    x->num_out_chans++;
  }
	
  x->num_in_chans = num_in_chans;
  alloc_memory(&x->temporary_memory, (t_uint) 0, (t_uint) 0);

  return (x);
}

t_convolve_error multi_channel_convolve_set(t_multi_channel_convolve *x, t_uint in_chan, t_uint out_chan, t_float *input, t_uint impulse_length, t_bool resize)
{
  /* For Parallel operation you must pass the same in/out channel */
	
  if(!x->N2M)
    in_chan -= out_chan;
	
  if(out_chan < x->num_out_chans)
    return output_channel_convolve_set(x->chan_convolvers[out_chan], in_chan, input, impulse_length, resize);
  else 
    return CONVOLVE_ERR_OUT_CHAN_OUT_OF_RANGE;
}

void multi_channel_convolve_clear(t_multi_channel_convolve *x, t_bool resize)
{
  t_uint i, j;
	
  if (x->N2M) {
    for (i = 0; i < x->num_out_chans; i++)
      for (j = 0; j < x->num_in_chans; j++)
	multi_channel_convolve_set(x, j, i, (t_float *)0, (t_uint) 0, resize);
  }
  else {
    for (i = 0; i < x->num_out_chans; i++)
      multi_channel_convolve_set(x, i, i, (t_float *)0, (t_uint) 0, resize);
  }
}


t_convolve_error multi_channel_convolve_resize(t_multi_channel_convolve *x, t_uint in_chan, t_uint out_chan, t_uint impulse_length)
{
  /* For Parallel operation you must pass the same in/out channel */
	
  if (!x->N2M)
    in_chan -= out_chan;
	
  if (out_chan < x->num_out_chans)
    return output_channel_convolve_resize(x->chan_convolvers[out_chan], in_chan, impulse_length);
  else 
    return CONVOLVE_ERR_IN_CHAN_OUT_OF_RANGE;
}


void multi_channel_convolve_process(t_multi_channel_convolve *x, t_sample **ins, t_sample **outs, t_sample *dry_gain, t_sample *wet_gain, t_uint vec_size, t_uint active_in_chans, t_uint active_out_chans)
{
  void *mem_pointer;
  vFloat **in_temps;
  vFloat *temp1;
  vFloat *temp2;
  t_uint num_in_chans = x->num_in_chans;
  t_uint num_out_chans = x->num_out_chans;
  long N2M = x->N2M;
  t_uint i, j;
	
#if defined( __i386__ ) || defined( __x86_64__ )
  unsigned int oldMXCSR = _mm_getcsr();
  /* read the old MXCSR setting */
  unsigned int newMXCSR = oldMXCSR | 0x8040;
  /* set DAZ and FZ bits */
  _mm_setcsr(newMXCSR);								  /* write the new MXCSR setting to the MXCSR */
#endif
	
  mem_pointer = grow_memory(&x->temporary_memory, (num_in_chans + 2) * vec_size * sizeof(t_float32), vec_size);  
  multi_channel_convolve_temp_setup(x, mem_pointer, x->temporary_memory.current_size);
	
  in_temps = x->in_temps;
  temp1 = x->temp1;
  temp2 = x->temp2;
	
  if(!x->temporary_memory.current_ptr)
    active_in_chans = active_out_chans = 0;
	
  active_in_chans = active_in_chans > num_in_chans ? num_in_chans : active_in_chans;
  active_out_chans = active_out_chans > num_out_chans ? num_out_chans : active_out_chans;
	
  for(i = 0; i < active_in_chans; i++) {
    t_float32 *current_in = (t_float32 *)in_temps[i];
		
    for (j = 0; j < vec_size; j++)
      current_in[j] = (t_float32)ins[i][j];
  }
	
  for (i = 0; i < active_out_chans; i++)
    output_channel_convolve_process(x->chan_convolvers[i], N2M ? in_temps : in_temps + i, outs[i], temp1, temp2, vec_size, active_in_chans);
	
  if(wet_gain) {
    for (i = 0; i < active_out_chans; i++)
      for (j = 0; j < vec_size; j++)
	outs[i][j] *= wet_gain[j];
  }
	
  if(dry_gain) {
    t_float32 *current_in;
		
    for (i = 0; i < active_in_chans && i < active_out_chans; i++)
      for (j = 0, current_in  = (t_float32 *)in_temps[i]; j < vec_size; j++)
	outs[i][j] += current_in[j] * dry_gain[j];
  }

	
#if defined( __i386__ ) || defined( __x86_64__ )	
  _mm_setcsr(oldMXCSR);	
#endif
}


