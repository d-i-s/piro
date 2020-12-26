#include "u_zero_latency_convolve.h"


void clear_memory(t_memory *mem_struct)
{
  if (mem_struct->current_free_method)
    mem_struct->current_free_method(mem_struct->current_ptr,
				    mem_struct->current_size);
}


void *equal_memory(t_memory *mem_struct, size_t size, size_t nom_size)
{	
  void *return_ptr;
  void *ptr;
		
  if (mem_struct->current_size != nom_size) {
    clear_memory(mem_struct);
			
    mem_struct->current_ptr = return_ptr = aligned_getbytes(size);
    mem_struct->current_size = return_ptr ? nom_size: 0;
    mem_struct->current_free_method = aligned_freebytes;
  }
  else
    return_ptr = mem_struct->current_ptr;
  
  return return_ptr;
}

void *equal_memory_custom(t_memory *mem_struct, alloc_method alloc_method_ptr, free_method free_method_ptr, size_t size, size_t nom_size)
{	
  void *return_ptr; 
	
  if (mem_struct->current_size != nom_size) {
    clear_memory(mem_struct);
		
    mem_struct->current_ptr = return_ptr = alloc_method_ptr(size, nom_size);
    mem_struct->current_size = return_ptr ? nom_size: 0;
    mem_struct->current_free_method = free_method_ptr;
  }
  else 
    return_ptr = mem_struct->current_ptr;
  
  return return_ptr;
}

void *attempt_memory_swap(t_memory *mem_struct, size_t *nom_size)
{
  void *return_ptr = NULL;
		
  *nom_size = mem_struct->current_size;
  return_ptr = mem_struct->current_ptr;
	
  return return_ptr;
}

void *grow_memory(t_memory *mem_struct, size_t size, size_t nom_size)
{	
  void *return_ptr;
	

  if (mem_struct->current_size < nom_size) {
    clear_memory(mem_struct);
		
    mem_struct->current_ptr = return_ptr = aligned_getbytes(size);
    mem_struct->current_size = return_ptr ? nom_size: 0;
    mem_struct->current_free_method = aligned_freebytes;
  }
  else 
    return_ptr = mem_struct->current_ptr;

  return return_ptr;
}

long alloc_memory(t_memory *mem_struct, size_t size, size_t nom_size)
{
  long fail = 0;
	
  if (size)
    mem_struct->current_ptr = aligned_getbytes(size);
  else
    mem_struct->current_ptr = NULL;

  if (size && mem_struct->current_ptr) {
    mem_struct->current_size = nom_size;
    mem_struct->current_free_method = aligned_freebytes;
  }
  else {
    mem_struct->current_size = 0;
    mem_struct->current_free_method = NULL;
		
    if (size) 
      fail = 1;
  }
	
  return fail;
}

long alloc_memory_custom(t_memory *mem_struct, alloc_method alloc_method_ptr, free_method free_method_ptr, size_t size, size_t nom_size)
{
  long fail = 0;
  
  if (size)
    mem_struct->current_ptr = alloc_method_ptr(size, nom_size);
  else
    mem_struct->current_ptr = NULL;
	
  if (size && mem_struct->current_ptr) {
    mem_struct->current_size = nom_size;
    mem_struct->current_free_method = free_method_ptr;
  }
  else {
    mem_struct->current_size = 0;
    mem_struct->current_free_method = NULL;
      
    if (size) 
      fail = 1;
  }
	
  return fail;
}

void *largest_partition_time_alloc(t_uint size, t_uint nom_size)
{
  return partition_convolve_new(16384, (size > 16384 ? size : 16384) - 8192, 8192, 0);
}


void *largest_partition_fft1_alloc(t_uint size, t_uint nom_size)
{
  return partition_convolve_new(16384, (size > 16384 ? size : 16384) - 8064, 8064, 0);
}

void *largest_partition_fft2_alloc(t_uint size, t_uint nom_size)
{
  return partition_convolve_new(16384, (size > 16384 ? size : 16384) - 7680, 7680, 0);
}


void largest_partition_free(void *large_part, size_t dummy)
{
  partition_convolve_free(large_part);
}


void zero_latency_convolve_free(t_zero_latency_convolve *x)
{
  if (!x)
    return;
	
  time_domain_convolve_free(x->time1);
  partition_convolve_free(x->part1);
  partition_convolve_free(x->part2);
  partition_convolve_free(x->part3);
  clear_memory(&x->part4);
  freebytes(x, sizeof(t_zero_latency_convolve));
}


t_zero_latency_convolve *zero_latency_convolve_new(t_uint max_length, t_convolve_latency_mode latency_mode)
{
  t_zero_latency_convolve *x = (t_zero_latency_convolve *)getbytes(sizeof(t_zero_latency_convolve));
  long fail = 0;
	
  if (!x)
    return 0;
	
  latency_mode = latency_mode > 2 ? 2 : latency_mode;

  size_t max_length_t = (size_t)max_length;

  switch (latency_mode) 
    {
    case CONVOLVE_LATENCY_ZERO:
		  
      x->time1 = time_domain_convolve_new(0, 128);
      x->part1 = partition_convolve_new(256, 384, 128, 384);
      x->part2 = partition_convolve_new(1024, 1536, 512, 1536);
      x->part3 = partition_convolve_new(4096, 6144, 2048, 6144);
      
      fail = alloc_memory_custom(&x->part4, largest_partition_time_alloc, largest_partition_free, max_length_t, max_length_t);
      
      break;
		
    case CONVOLVE_LATENCY_SHORT:
			
      x->time1 = NULL;
      x->part1 = partition_convolve_new(256, 384, 0, 384);
      x->part2 = partition_convolve_new(1024, 1536, 384, 1536);
      x->part3 = partition_convolve_new(4096, 6144, 1920, 6144);
      fail = alloc_memory_custom(&x->part4, largest_partition_fft1_alloc, largest_partition_free, max_length_t, max_length_t);
			
      break;
			
    case CONVOLVE_LATENCY_MEDIUM:
			
      x->time1 = NULL;
      x->part1 = NULL;
      x->part2 = partition_convolve_new(1024, 1536, 0, 1536);
      x->part3 = partition_convolve_new(4096, 6144, 1536, 6144);
      fail = alloc_memory_custom(&x->part4, largest_partition_fft2_alloc, largest_partition_free, max_length_t, max_length_t);
			
      break;
			
    }

  x->latency_mode = latency_mode;
  x->impulse_length = 0;
	
  if ((!latency_mode && !x->time1) || (latency_mode < 2 && !x->part1) || !x->part2 || !x->part3 || fail) {
    zero_latency_convolve_free(x);
    return 0;
  }
	
  return (x);
}


t_partition_convolve *zero_latency_convolve_resize(t_zero_latency_convolve *x, t_uint impulse_length, t_bool keep_lock)
{
  t_partition_convolve *return_part = NULL;
  x->impulse_length = 0;
	
  switch (x->latency_mode)
    {
    case 0:
      return_part = equal_memory_custom(&x->part4, largest_partition_time_alloc, largest_partition_free, impulse_length, impulse_length);	
      break;
			
    case 1:
      return_part = equal_memory_custom(&x->part4, largest_partition_fft1_alloc, largest_partition_free, impulse_length, impulse_length);	
      break;

    case 2:
      return_part = equal_memory_custom(&x->part4, largest_partition_fft2_alloc, largest_partition_free, impulse_length, impulse_length);	
      break;
    }	
	
	
  return return_part;
}


t_convolve_error zero_latency_convolve_set(t_zero_latency_convolve *x, t_float *input, t_uint impulse_length, t_bool resize)
{	
  t_partition_convolve *part4 = NULL;
  t_uint max_impulse = 0;
	
  x->impulse_length = 0;
  
  if(resize) {
    part4 = zero_latency_convolve_resize(x, impulse_length, true);
    max_impulse = impulse_length;
  }
  else 
    part4 = (t_partition_convolve *)x->part4.current_ptr;

  if(part4) {
    if (x->latency_mode < 1)
      time_domain_convolve_set(x->time1, input, impulse_length);
    if (x->latency_mode < 2)
      partition_convolve_set(x->part1, input, impulse_length);
    partition_convolve_set(x->part2, input, impulse_length);
    partition_convolve_set(x->part3, input, impulse_length);
    partition_convolve_set(part4, input, impulse_length);
  }

  x->impulse_length = impulse_length;
	
  if(impulse_length && !part4)
    return CONVOLVE_ERR_MEM_UNAVAILABLE;
	
  if(impulse_length > max_impulse)
    return CONVOLVE_ERR_MEM_ALLOC_TOO_SMALL;
		
  return CONVOLVE_ERR_NONE;
}


void zero_latency_convolve_process_sum(vFloat *out, vFloat *add, t_uint vec_size)
{
  t_uint i;
	
  for (i = 0; i < (vec_size >> 2); i++, out++) 
    *out = F32_VEC_ADD_OP(*out, *add++);
  
}


void zero_latency_convolve_process(t_zero_latency_convolve *x, vFloat *in, vFloat *temp, vFloat *out, t_uint vec_size)
{		
  size_t max_impulse = 0; 
  t_partition_convolve *part4 = attempt_memory_swap(&x->part4, &max_impulse);
  
  /* N.B. This function DOES NOT zero the output buffer as this is done elsewhere */
	
  if (x->impulse_length && x->impulse_length <= max_impulse && part4) {

    if (x->latency_mode == 0) {
      time_domain_convolve_process(x->time1, (t_float32 *)in, (t_float32 *)temp,
				   vec_size);
      zero_latency_convolve_process_sum(out, temp, vec_size);
      
    }
    if (x->latency_mode < 2) {
      if (partition_convolve_process(x->part1, in, temp, vec_size) == true)
	zero_latency_convolve_process_sum(out, temp, vec_size);
    }
    if (partition_convolve_process(x->part2, in, temp, vec_size) == true)
      zero_latency_convolve_process_sum(out, temp, vec_size);
    if (partition_convolve_process(x->part3, in, temp, vec_size) == true)
      zero_latency_convolve_process_sum(out, temp, vec_size);
    if (partition_convolve_process(part4, in, temp, vec_size) == true)
      zero_latency_convolve_process_sum(out, temp, vec_size);
  }
	
}

