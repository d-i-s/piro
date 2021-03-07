#ifndef __Z_CORE__
#define __Z_CORE__

#include "z_fft.h"

#define PIRO_MAX_SPECIFIER_ITEMS 2048
#define PIRO_DB_MIN -500

#if !defined(HUGE_VAL)
#error HUGE_VAL not defined
#define HUGE_VAL 7FFFFFFF /* ridefiniamo questo valore per i 32 bit? */
#endif

static t_float hann_table[4097]; /* thanks to LC */
static long hann_setup_flag = 0;

typedef enum {
  SPECTRUM_REAL = 0,
  SPECTRUM_FULL = 1
} t_spectrum_format;

typedef enum {
  FILTER_REGULARISATION = 0,
  FILTER_CLIP = 1,
  FILTER_FILTER = 2
} t_filter_type;

double db_to_pow(double db);
void db_to_pow_array(t_sample *in, t_uint length);
double pow_to_db(double pow);
double db_to_a(double db);


void spike_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format, double spike);
void delay_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format, double delay);


void make_freq_dependent_power_array(t_sample *power_array, t_sample *specifier_array, t_uint fft_size, t_sample sample_rate, t_sample db_offset);


void deconvolve_with_amp_filter(FFT_Split spectrum_1, FFT_Split spectrum_2, t_float *filter_amps, t_uint fft_size, t_spectrum_format format);
void deconvolve_clip_zero_phase(FFT_Split spectrum_1, FFT_Split spectrum_2, t_float *clip_min, t_float *clip_max, t_uint fft_size, t_spectrum_format format);
void deconvolve_regularised_zero_phase(FFT_Split spectrum_1, FFT_Split spectrum_2, t_float *beta_in, t_uint fft_size, t_spectrum_format format);
void deconvolve_zero_phase(FFT_Split spectrum_1, FFT_Split spectrum_2, FFT_Split filter_spectrum, t_float *filter_specifier, t_float *range_specifier, t_float filter_db_offset, t_uint fft_size, t_spectrum_format format, t_filter_type mode, t_float sample_rate);
void deconvolve_with_filter(FFT_Split spectrum_1, FFT_Split spectrum_2, FFT_Split filter_spectrum, t_uint fft_size, t_spectrum_format format);


void zero_phase_from_power_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format);
void linear_phase_from_power_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format);
void minimum_phase_components_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size);
void minimum_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size);
void noncausal_maximum_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size);
void maximum_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size);
void mixed_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size, double phase, t_bool zero_center);
void variable_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size, double phase, t_bool zero_center);


void power_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format);
void make_regularisation_filter(FFT_Setup *fft_setup, FFT_Split denominator_spectrum, FFT_Split filter_spectrum, t_sample *beta_in, t_uint fft_size, t_spectrum_format format, double phase);
void make_clip_filter(FFT_Setup *fft_setup, FFT_Split denominator_spectrum, FFT_Split filter_spectrum, t_sample *clip_min, t_sample *clip_max, t_uint fft_size, t_spectrum_format format, double phase);
void time_to_spectrum(FFT_Setup *fft_setup, t_sample *in_buf, t_uint in_length, FFT_Split spectrum, t_uint fft_size);
void make_deconvolution_filter(FFT_Setup *fft_setup, FFT_Split denominator_spectrum, FFT_Split filter_spectrum, t_float *filter_specifier, t_float *range_specifier, t_sample filter_db_offset, t_float *filter_in, t_uint filter_length, t_uint fft_size, t_spectrum_format format, t_filter_type mode, double phase, double sample_rate);
void deconvolve_variable_phase(FFT_Setup *fft_setup, FFT_Split spectrum_1, FFT_Split spectrum_2, FFT_Split filter_spectrum, t_sample *filter_specifier, t_sample *range_specifier, t_float filter_db_offset, t_float *filter_in, t_uint filter_length, t_uint fft_size, t_spectrum_format format, t_filter_type mode, t_float phase, t_float sample_rate);
void deconvolve(FFT_Setup *fft_setup, FFT_Split spectrum_1, FFT_Split spectrum_2, FFT_Split filter_spectrum, t_float *filter_specifier, t_float *range_specifier, t_sample filter_db_offset, t_float *filter_in, t_uint filter_length, t_uint fft_size, t_spectrum_format format, t_filter_type mode, t_float phase, t_float delay, t_float sample_rate);
void time_to_halfspectrum(FFT_Setup *fft_setup, t_sample *in_buf, t_uint in_length, FFT_Split spectrum, t_uint fft_size);
void convolve(FFT_Split fft_data_1, FFT_Split fft_data_2, t_uint fft_size,
	      t_spectrum_format format);
void power_full_spectrum_from_half_spectrum(FFT_Split spectrum, t_uint fft_size);
void spectrum_to_time(FFT_Setup *fft_setup, t_sample *out_buf, FFT_Split spectrum, t_uint fft_size, t_spectrum_format half_spectrum);


void setup_hann_wind();
t_float fast_hann_wind(t_float in);
void smooth_power_spectrum(FFT_Split spectrum, long mode, t_uint fft_size, t_float smooth_lo, t_float smooth_hi);


/* #ifndef __APPLE__ */

/* static double round(double r) */
/* { */
/*   return (r > 0.) ? floor(r + 0.5) : ceil(r - 0.5); */
/* } */

/* static long isnan(double n) */
/* { */
/*   return !(n == n); */
/* } */

/* static long isinf(double n) */
/* { */
/*   return !isnan(n) & isnan(n - n); */
/* } */

/* #endif */ /* APPLE */

#endif /* __ZCORE__ */
