#ifndef __Z_FFT__
#define __Z_FFT__


#include "m_types.h"

#if defined (NT)
#include <intrin.h>
#endif

static int SSE_Exists = 0;

/* static int SSE2_check(); */

void pass_1_2_reorder_simd32(FFT_Split32 *input, t_uint length);
void pass_3_reorder_simd32(FFT_Split32 *input, t_uint length);
void pass_trig_table_reorder_simd32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint length, t_uint pass);
void pass_trig_table_simd32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint length, t_uint pass);

void pass_real_trig_table32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2, long ifft);
void do_small_real_fft32(FFT_Split32 *input, t_uint fft_log2, long ifft);
void pass_trig_table32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint length, t_uint pass);
void pass_trig_table_reorder32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint length, t_uint pass);
void pass_3_reorder32(FFT_Split32 *input, t_uint length, t_uint fft_log2);
void pass_1_2_reorder32(FFT_Split32 *input, t_uint length, t_uint fft_log2);
void pass_3_32(FFT_Split32 *input, t_uint length, t_uint fft_log2);
void do_small_fft32(FFT_Split32 *input, t_uint fft_log2);


void pass_1_2_reorder_simd(FFT_Split *input, t_uint length);
void pass_3_reorder_simd(FFT_Split *input, t_uint length);
void pass_trig_table_reorder_simd(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass);
void pass_trig_table_simd(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass);

void pass_real_trig_table(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2, long ifft);
void do_small_real_fft(FFT_Split *input, t_uint fft_log2, long ifft);
void pass_trig_table(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass);
void pass_trig_table_reorder(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass);
void pass_3_reorder(FFT_Split *input, t_uint length, t_uint fft_log2);
void pass_1_2_reorder(FFT_Split *input, t_uint length, t_uint fft_log2);
void pass_3(FFT_Split *input, t_uint length, t_uint fft_log2);
void do_small_fft(FFT_Split *input, t_uint fft_log2);

t_uint int_log2(t_uint in, t_uint *inexact);
t_uint calculate_fft_size(t_uint input_size, t_uint *fft_size_log2);


void fft_fill_table(FFT_Split *table, t_uint length);
FFT_Setup *create_setup(t_uint max_fft_log2);
void destroy_setup(FFT_Setup *setup);
void fft_fill_table32(FFT_Split32 *table, t_uint length);
FFT_Setup32 *create_setup32(t_uint max_fft_log2);
void destroy_setup32(FFT_Setup32 *setup);


void fft_internal32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2);
void do_fft32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2);
void do_ifft32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2);
void fft_internal(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2);
void do_fft(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2);
void do_ifft(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2);

void do_real_fft(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2);
void do_real_ifft(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2);
void do_real_fft32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2);
void do_real_ifft32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2);

void unzip_complex32(t_float32 *input, FFT_Split32 *output, t_uint half_length);
void unzip_complex(t_sample *input, FFT_Split *output, t_uint half_length);
void unzip_zero(t_sample *input, FFT_Split *output, t_uint in_length, t_uint log2n);
void zip_sample32(FFT_Split32 *input, t_float32 *output, t_uint half_length);
void zip_sample(FFT_Split *input, t_sample *output, t_uint half_length);


#if defined( __x86_64__ ) || defined( __LP64__ )  || defined(_WIN64)
typedef unsigned long long t_uint_fft;
#else
typedef unsigned long t_uint_fft;
#endif

#endif /* Z_FFT */
