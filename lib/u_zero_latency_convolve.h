#include "m_types.h"

#include "u_partition_convolve.h"
#include "u_time_domain_convolve.h"
#include "m_convolve_latency_modes.h"
#include "m_convolve_errors.h"

#ifndef __U_ZEROLATENCYCONVOLVE_STRUCT__
#define __U_ZEROLATENCYCONVOLVE_STRUCT__


typedef void *(*alloc_method) (t_uint size, t_uint nom_size);
typedef void (*free_method) (void *, size_t size);



/* free method memory structure */

typedef struct _memory
{
  void *current_ptr;
  free_method current_free_method;
  size_t current_size;	
} t_memory;

void clear_memory(t_memory *mem_struct);
void *equal_memory(t_memory *mem_struct, size_t size, size_t nom_size);
void *equal_memory_custom(t_memory *mem_struct, alloc_method alloc_method_ptr, free_method free_method_ptr, size_t size, size_t nom_size);
void *grow_memory(t_memory *mem_struct, size_t size, size_t nom_size);
long alloc_memory(t_memory *mem_struct, size_t size, size_t nom_size);
long alloc_memory_custom(t_memory *mem_struct, alloc_method alloc_method_ptr, free_method free_method_ptr, size_t size, size_t nom_size);



typedef struct _zero_latency_convolve
{
  t_time_domain_convolve *time1;
  t_partition_convolve *part1;
  t_partition_convolve *part2;
  t_partition_convolve *part3;
	
  t_memory part4;
	
  t_uint impulse_length;
  t_convolve_latency_mode latency_mode;
	
} t_zero_latency_convolve;

#endif // __U_ZEROLATENCYCONVOLVE_STRUCT__

void zero_latency_convolve_free(t_zero_latency_convolve *x);
t_zero_latency_convolve *zero_latency_convolve_new(t_uint max_length, t_convolve_latency_mode latency_mode);

t_partition_convolve *zero_latency_convolve_resize(t_zero_latency_convolve *x, t_uint impulse_length, t_bool keep_lock);
t_convolve_error zero_latency_convolve_set(t_zero_latency_convolve *x, t_float *input, t_uint impulse_length, t_bool resize);

void zero_latency_convolve_process_sum(vFloat *out, vFloat *add, t_uint vec_size);
void zero_latency_convolve_process(t_zero_latency_convolve *x, vFloat *in, vFloat *temp, vFloat *out, t_uint vec_size);

