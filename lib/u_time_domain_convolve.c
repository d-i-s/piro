
/*
 *  time_domain_convolve
 *
 *	time_domain_convolve performs real-time zero latency time-based convolution.
 *	
 *	Typically time_domain_convolve is suitable for use in conjunction with partition_convolve for zero-latency convolution with longer impulses (time_domain_convolve use apple's vDSP and the IR length is limited to 2044 samples).
 *	Note that in fact the algorithms process correlation with reversed impulse response coeffients - which is equivalent to convolution.
 *
 *  Copyright 2012 Alex Harker. All rights reserved.
 *
 */


#include "u_time_domain_convolve.h"


t_uint pad_length(t_uint length)
{
  return ((length + 15) >> 4) << 4;
}


void time_domain_convolve_offset_set(t_time_domain_convolve *x, t_uint offset)
{
  x->offset = offset;
}


t_convolve_error time_domain_convolve_length_set(t_time_domain_convolve *x, t_uint length)
{	
  t_convolve_error error = CONVOLVE_ERR_NONE;

  if (length > 2044) {
    error = CONVOLVE_ERR_TIME_LENGTH_OUT_OF_RANGE;
    length = 2044;
  }
	
  x->length = length;
	
  return error;
}

void time_domain_convolve_free(t_time_domain_convolve *x)
{
  if (!x)
    return;
	
  aligned_freebytes(x->impulse_buffer, sizeof(t_float32)*2048);
  aligned_freebytes(x->input_buffer, sizeof(t_float32)*8192);
  freebytes(x, sizeof(t_time_domain_convolve));
}

t_time_domain_convolve *time_domain_convolve_new(t_uint offset, t_uint length)
{
  t_time_domain_convolve *x = (t_time_domain_convolve *)getbytes(sizeof(t_time_domain_convolve));
  t_uint i;

  if (!x)
    return 0;
	
  /* Set default initial variables */
	
  x->input_position = 0;
  x->impulse_length = 0;
	
  time_domain_convolve_offset_set(x, offset);
  time_domain_convolve_length_set(x, length);
	
  /* Allocate impulse buffer and input buffer */
	
  x->impulse_buffer = aligned_getbytes(sizeof(t_float32) * 2048);
  x->input_buffer = aligned_getbytes(sizeof(t_float32) *  8192);
	
  if (!x->impulse_buffer || !x->input_buffer) {
    time_domain_convolve_free(x);
    return 0;
  }
	
  for (i = 0; i < 2048; i++) /* non dovrebbe servire */
    x->impulse_buffer[i] = 0.f;

  for (i = 0; i < 8192; i++)
    x->input_buffer[i] = 0.f;

	
  return (x);
}


t_convolve_error time_domain_convolve_set(t_time_domain_convolve *x, t_float *input, t_uint impulse_length)
{	
  t_uint offset = x->offset;
  t_uint length = x->length;
	
  t_float32 *impulse_buffer = x->impulse_buffer;

  t_convolve_error error = CONVOLVE_ERR_NONE;

#if defined (NT)
  t_uint impulse_offset;
#endif
  
  t_uint i, j;
	
  x->impulse_length = 0;
	
  /* Calculate impulse length */
	
  if (!input || impulse_length < offset)
    impulse_length = 0;
	
  impulse_length -= offset;
  if (length && length < impulse_length)
    impulse_length = length;
	
  if (impulse_length > 2044) {
    error = CONVOLVE_ERR_TIME_IMPULSE_TOO_LONG;
    impulse_length = 2044;
  }

#if defined (UNIX) || defined (MACOSX)
  if (impulse_length) {
    for (i = impulse_length, j = 0; i > 0; i--, j++)
      impulse_buffer[j] = (t_float32)input[i + offset - 1];
  }
#else
  if (impulse_length) {
    impulse_offset = pad_length(impulse_length) - impulse_length; 

    for (i = 0; i < impulse_offset; i++)
      impulse_buffer[i] = 0.f;

    for (i = impulse_length, j = 0; i > 0; i--, j++)
      impulse_buffer[j + impulse_offset] = (t_float32)input[i + offset - 1];
  }
#endif
		
  x->impulse_length = impulse_length;
	
  return error;
}

void td_conv(t_float32 *in, vFloat *impulse, t_float32 *output, t_uint N, t_uint L)
{
  vFloat output_accum;
  t_float32 *input;
  t_float32 results[4];

  t_uint i, j;
		
  L = pad_length(L);
				   
  for (i = 0; i < N; i++) {
    output_accum = float2vector(0.f);
    input = in - L + 1 + i;

    for (j = 0; j < L >> 2; j += 4) {
      /* Load vals */
			
      output_accum = F32_VEC_ADD_OP(output_accum, F32_VEC_MUL_OP(impulse[j], F32_VEC_ULOAD(input)));
      input += 4;
      output_accum = F32_VEC_ADD_OP(output_accum, F32_VEC_MUL_OP(impulse[j + 1], F32_VEC_ULOAD(input)));
      input += 4;
      output_accum = F32_VEC_ADD_OP(output_accum, F32_VEC_MUL_OP(impulse[j + 2], F32_VEC_ULOAD(input)));
      input += 4;
      output_accum = F32_VEC_ADD_OP(output_accum, F32_VEC_MUL_OP(impulse[j + 3], F32_VEC_ULOAD(input)));
      input += 4;
    }

    F32_VEC_USTORE(results, output_accum);
    *output++ = results[0] + results[1] + results[2] + results[3];
  }
}


void td_conv_scalar(t_float32 *in, t_float32 *impulse, t_float32 *output, t_uint N, t_uint L)
{
  t_float32 output_accum;
  t_float32 *input;
	
  t_uint i, j;
	
  L = pad_length(L);
	
  for (i = 0; i < N; i++) {
    output_accum = 0.f;
    input = in - L + 1 + i;
		
    for (j = 0; j < L; j += 8) {
      /* Load vals */
			
      output_accum += impulse[j+0] * *input++;
      output_accum += impulse[j+1] * *input++;
      output_accum += impulse[j+2] * *input++;
      output_accum += impulse[j+3] * *input++;
      output_accum += impulse[j+4] * *input++;
      output_accum += impulse[j+5] * *input++;
      output_accum += impulse[j+6] * *input++;
      output_accum += impulse[j+7] * *input++;
    }
		
    *output++ = output_accum;
  }
}


void time_domain_convolve_process_scalar(t_time_domain_convolve *x, t_float32 *in, t_float32 *out, t_uint vec_size)
{
  t_float32 *impulse_buffer = x->impulse_buffer;
  t_float32 *input_buffer = x->input_buffer;
  t_uint input_position = x->input_position;
  t_uint impulse_length = x->impulse_length;
  t_uint current_loop;

  for (current_loop = vec_size > 4096 ? 4096 : vec_size; current_loop; vec_size -= current_loop, current_loop = vec_size > 4096 ? 4096 : vec_size) {
    
  /* Copy input twice (allows us to read input out in one go) */
	
  memcpy(input_buffer + input_position, in, sizeof(t_float32) * vec_size);
  memcpy(input_buffer + 4096 + input_position, in, sizeof(t_float32) * vec_size);
	
  /* Advance pointer */
	
  input_position += vec_size;
  if (input_position >= 4096) 
    input_position -= 4096;
  x->input_position = input_position;

  t_float32 *tmp = input_buffer + 4096 + input_position - vec_size;
	
  /* Do convolution */
	
  td_conv_scalar(tmp, impulse_buffer, out, vec_size, impulse_length);
  
  }
}



void time_domain_convolve_process(t_time_domain_convolve *x, t_float32 *in, t_float32 *out, t_uint vec_size)
{	
  t_float32 *impulse_buffer = x->impulse_buffer;
  t_float32 *input_buffer = x->input_buffer;
  t_uint input_position = x->input_position;
  t_uint impulse_length = x->impulse_length;
  t_uint current_loop;
  
  for (current_loop = vec_size > 4096 ? 4096 : vec_size; current_loop; vec_size -= current_loop, current_loop = vec_size > 4096 ? 4096 : vec_size) {
    /* Copy input twice (allows us to read input out in one go) */
		
    memcpy(input_buffer + input_position, in, sizeof(t_float32) * vec_size);
    memcpy(input_buffer + 4096 + input_position, in, sizeof(t_float32) * vec_size);
		
    /* Advance pointer */
		
    input_position += vec_size;
    if (input_position >= 4096) 
      input_position -= 4096;
    x->input_position = input_position;
		
    /* Do convolution */

    t_float32 *tmp = input_buffer + 4096 + input_position - vec_size;

    td_conv(tmp, (vFloat *) impulse_buffer, out, vec_size, impulse_length);

  }
}

