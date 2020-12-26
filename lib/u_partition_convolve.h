#include "m_convolve_errors.h"

#include "z_fft.h"
#include "m_types.h"

#ifdef __linux__
#define arc4random rand
#endif /* __linux__ */

#ifdef NT
#include <Windows.h>
#endif

#ifndef __U_PARTCONVOLVE_STRUCT__
#define __U_PARTCONVOLVE_STRUCT__

typedef struct _partition_convolve
{		
  /* FFT variables */
	
  FFT_Setup32 *fft_setup_real;
	
  t_uint max_fft_size;
  t_uint max_fft_size_log2; 
  t_uint fft_size; 
  t_uint fft_size_log2;
	
  t_uint till_next_fft; 
  t_uint rw_pointer1;
  t_uint rw_pointer2;
	
  /* Scheduling variables */
	
  t_uint num_partitions;
  t_uint valid_partitions;
  t_uint partitions_done;
  t_uint last_partition;
	
  t_uint input_position;
  t_uint schedule_counter;
	
  /* Internal buffers */
	
  vFloat *fft_buffers[4];
	
  FFT_Split32 impulse_buffer;
  FFT_Split32 input_buffer;
  FFT_Split32 accum_buffer;
  FFT_Split32 partition_temp;
	
  t_uint max_impulse_length;
	
  /* Attributes */
	
  t_uint offset;
  t_uint length;
	
  /* Flags */
	
  char reset_flag; /* reset fft data on next perform call */
	
} t_partition_convolve;

#endif /* __U_PARTCONVOLVE_STRUCT__ */


// N.B. MIN_FFT_SIZE_LOG2 should never be smaller than 4, as below code assumes loop unroll of vectors (4 vals) by 4 (== 16 or 2^4)
//		MAX_FFT_SIZE_LOG2 is perhaps conservative right now, assuming realtime usage, but it is easy to increase this if necessary

#define MIN_FFT_SIZE_LOG2					5
#define MAX_FFT_SIZE_LOG2					20


void partition_convolve_free(t_partition_convolve *x);
t_partition_convolve *partition_convolve_new(t_uint max_fft_size, t_uint max_impulse_length, t_uint offset, t_uint length);
void init_partition_convolve();

t_convolve_error partition_convolve_fft_size_set(t_partition_convolve *x, t_uint fft_size);
t_convolve_error partition_convolve_length_set(t_partition_convolve *x, t_uint length);
void partition_convolve_offset_set(t_partition_convolve *x, t_uint offset);

t_convolve_error partition_convolve_set(t_partition_convolve *x, t_float *input, t_uint impulse_length);

t_bool partition_convolve_process(t_partition_convolve *x, vFloat *in, vFloat *out, t_uint vec_size);
void partition_convolve_process_partition(FFT_Split32 in1, FFT_Split32 in2, FFT_Split32 out, t_uint num_vecs);
