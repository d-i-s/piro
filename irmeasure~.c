#include "irmeasure~.h"

typedef enum {
  SWEEP = 0,
  MLS = 1,
  NOISE = 2
} t_excitation_signal;

static t_class *irmeasure_tilde_class;

typedef struct _irmeasure_tilde {
  t_object  x_obj;
  t_sample x_f;

  /* Requested Sweep / MLS / Noise Parameters */
	
  double lo_f;
  double hi_f;
  double fade_in;
  double fade_out;
  double length;
  double out_length;
  double sample_rate;
		
  long order;
  long num_active_ins;
  long num_active_outs;
	
  t_noise_mode noise_mode;
	
  /* Internal */
	
  double test_tone_freq;
  double current_out_length;
  double phase;
	
  long start_measurement;
  long stop_measurement;
  long test_tone;
  long no_dsp;
	
  t_int current_t;
  t_int T;
  t_int T2;
  t_int fft_size;
  t_int chan_offset[PIRO_MAX_MEASURE_CHANS]; // PIRO_
	
  long num_in_chans;
  long num_out_chans;
  long current_num_active_ins;
  long current_num_active_outs;

  void *in_chans[PIRO_MAX_MEASURE_CHANS]; // PIRO_
  void *out_chans[PIRO_MAX_MEASURE_CHANS + 1]; // PIRO_
	
  /* Noise Amp Compensation */
	
  double max_amp_brown;
  double max_amp_pink;
	
  /* Measurement Parameters */
	
  t_excitation_signal measure_mode;
	
  t_ess sweep_params;
  t_mls max_length_params;
  t_noise_params noise_params;
	
  /* Permanent Memory */
	
  t_sample *rec_mem;
  t_sample *out_mem;
  t_int rec_mem_size;
  t_int out_mem_size;
	
  /* Amplitude Curve Temp */
	
  double amp_curve[33];

  /* attributes */
  double amp; /* signal amplitude (dB) */
  long abs_progress; /* absolute progress */
  long bandlimit; /* bandlimit sweep measurements */
  long inv_amp; /* invert amplitude */

  /* sweep attr */
  t_atom amp_curve_specifier[32];
  long amp_curve_num_specifiers;

  /* deconvolution attr */
  t_atom *deconvolve_filter_specifier;
  t_atom *deconvolve_range_specifier;
  long deconvolve_num_filter_specifiers;
  long deconvolve_num_range_specifiers;
  long deconvolve_mode;
  t_atom deconvolve_phase;
  t_atom deconvolve_delay; // non necessario

  /* buffer write attr */
  long resize;
  long write_chan;
  
  t_outlet *bangout; // free
  
} t_irmeasure_tilde;

static const t_uint32 feedback_mask_vals[] = {0x0u, 0x0u, 0x2u, 0x6u, 0xCu, 0x14u, 0x30u, 0x60u, 0xE1u, 0x100u, 0x240u, 0x500u, 0xE08u, 0x1C80u, 0x3802u, 0x6000u, 0xD008, 0x12000u, 0x20400u, 0x72000u, 0x90000u, 0x500000u, 0xC00000u, 0x420000u, 0xE10000u};

void fill_amp_curve_specifier(double *array, t_atom *specifier, long num_specifiers)
{
  long i=0;

  for(i=0;i<num_specifiers;i++) {
    array[i] = atom_getfloat(specifier+i);
  }

  array[i] = HUGE_VAL;
}

void fill_power_array_specifier(t_float *array, t_atom *specifier, long num_specifiers)
{
  long i;

  num_specifiers = num_specifiers > PIRO_MAX_SPECIFIER_ITEMS ?
    PIRO_MAX_SPECIFIER_ITEMS : num_specifiers;

  if(specifier->a_type == A_SYMBOL) {
    array[0] = -100.;
    array[1] = HUGE_VAL;
    return;
  }
  
  for(i=0; i<num_specifiers; i++)
    array[i] = atom_getfloat(specifier+i);

  if(num_specifiers < PIRO_MAX_SPECIFIER_ITEMS)
    array[i] = HUGE_VAL;
}

static t_uint coloured_noise_get_length(t_noise_params *x)
{
  return x->T;
}

static t_uint mls_get_length(t_mls *x)
{
  return x->T;
}

static t_uint ess_get_length(t_ess *x)
{
  return x->T;
}

t_int irmeasure_calc_sweep_mem_size(t_ess *sweep_params, long num_out_chans, double out_length, double sample_rate)
{
  t_int gen_length = ess_get_length(sweep_params);
  t_int rec_length = (t_int) (num_out_chans * ((out_length * sample_rate) + gen_length));

  return rec_length * sizeof(t_sample);
}

t_int irmeasure_calc_mls_mem_size(long order, long num_out_chans, double out_length, double sample_rate)
{
  t_int gen_length = (1 << order) - 1;
  t_int rec_length = (t_int) (num_out_chans * ((out_length *
							  sample_rate) +
							 gen_length));

  return rec_length * sizeof(t_sample);
}

t_int irmeasure_calc_noise_mem_size(double length, long num_out_chans, double out_length, double sample_rate)
{
  t_int gen_length = (t_int) (length * sample_rate);
  t_int rec_length = (t_int) (num_out_chans * ((out_length *
							  sample_rate) +
							 gen_length));

  return rec_length * sizeof(t_sample);
}

t_int irmeasure_calc_mem_size(t_irmeasure_tilde *x, long num_in_chans,
			   long num_out_chans, double sample_rate)
{
  switch(x->measure_mode)
    {
    case SWEEP:
      return num_in_chans * irmeasure_calc_sweep_mem_size(&x->sweep_params,
							  num_out_chans,
							  x->current_out_length,
							  sample_rate);
    case MLS:
      return num_in_chans * irmeasure_calc_mls_mem_size(x->max_length_params.order, num_out_chans, x->current_out_length, sample_rate);
    case NOISE:
      return num_in_chans * irmeasure_calc_noise_mem_size(x->length,
							  num_out_chans,
							  x->current_out_length,
							  sample_rate);
    }
  return 0;
}

void coloured_noise_reset(t_noise_params *x)
{
  x->gen.w = 0;
  x->gen.x = 0;
  x->gen.y = 0;
  x->gen.z = 4294967295u;
}

void coloured_noise_params(t_noise_params *x, t_noise_mode mode, double fade_in, double fade_out, double length, double sample_rate, double amp)
{
  x->alpha = sin(M_PI * 2.0 * 16.0 / sample_rate);

  x->alpha0 = sin(M_PI * 2.0 * 8.00135734209627 / sample_rate);
  x->alpha1 = sin(M_PI * 2.0 * 46.88548507044182 / sample_rate);
  x->alpha2 = sin(M_PI * 2.0 * 217.61558695916962 / sample_rate);
  x->alpha3 = sin(M_PI * 2.0 * 939.80665948455472 / sample_rate);
  x->alpha4 = sin(M_PI * 2.0 * 3276.10128392439381 / sample_rate);

  x->prev_output = 0.;

  x->b0 = 0.;
  x->b1 = 0.;
  x->b2 = 0.;
  x->b3 = 0.;
  x->b4 = 0.;
  x->b5 = 0.;
  x->b6 = 0.;

  coloured_noise_reset(x);

  x->amp = amp;
  x->sample_rate = sample_rate;
  x->T = (t_uint)(length * sample_rate);
  x->RT = length;
  x->fade_in = fade_in;
  x->fade_out = fade_out;

  if(mode>2)
    mode = 0;

  x->mode = mode;
}

void mls_reset(t_mls *x)
{
  x->lfsr = 0x1u;
}

void mls_params(t_mls *x, t_uint32 log2_T, double amp)
{
  log2_T = (log2_T < 1) ? 1 : log2_T;
  log2_T = (log2_T > 24) ? 24 : log2_T;
  x->amp = amp;

  x->T = (1u << log2_T) - 1u;
  x->order = log2_T;
  x->feedback_mask = feedback_mask_vals[log2_T];
  x->lfsr = 0x1u;
}

t_uint ess_params(t_ess *x, double f1, double f2, double fade_in,
		  double fade_out, double T, double sample_rate, double amp,
		  double *amp_curve)
{
  double L;

  double K1;
  double K2;

  double nt;
  double final_phase;
  double NNT;

  double last_db_val = 0.;

  unsigned long num_items = 0;
  unsigned int i = 0;

  x->RT = T;
  x->rf1 = f1;
  x->rf2 = f2;

  f1 /= sample_rate;
  f2 /= sample_rate;

  T *= sample_rate;

  L = round(f1 * T / (log(f2/f1))) / f1;

  K1 = 2 * M_PI * f1 * L;
  K2 = 1 / L;

  nt = round(f1 * T / (log(f2/f1))) * log(f2/f1) / f1;
  final_phase = floor(L * f1 * (exp(nt * K2) - 1));
  NNT = ceil(log((final_phase / L / f1 + 1)) / K2);

  if(L==0)
    return 0;

  x->K1 = K1;
  x->K2 = K2;
  x->T = (t_uint)NNT;
  x->lo_f_act = f1;
  x->hi_f_act = f1 * exp(NNT / L);
  x->f1 = f1;
  x->f2 = f2;
  x->fade_in = fade_in;
  x->fade_out = fade_out;
  x->sample_rate = sample_rate;
  x->amp = amp;
  
  for(i=0; amp_curve && (i < 32); i++) { /* isinf fix */
    if(isinf(amp_curve[i]) || amp_curve==NULL || amp_curve[i] == HUGE_VAL) 
      break;
  }
  
  if(i!=0) {
    num_items = i >> 1;
  }

  x->amp_specifier[0] = 0.f;
  x->amp_specifier[1] = num_items ? amp_curve[1] : 0.f;

  for(i=0; i<num_items; i++) {
    x->amp_specifier[2*i + 2] = log(amp_curve[2*i] / (f1 * sample_rate));
    x->amp_specifier[2*i + 3] = last_db_val = amp_curve[2*i + 1];
  }

  x->amp_specifier[2*num_items + 2] = HUGE_VAL;
  x->amp_specifier[2*num_items + 3] = last_db_val;

  x->num_amp_specifiers = num_items;
  
  return (t_uint)NNT;
}

void irmeasure_sweep_params(t_irmeasure_tilde *x)
{
  double out_length = x->out_length;
  double sample_rate = x->sample_rate;

  long sweep_length = (long)ess_params(&x->sweep_params, x->lo_f, x->hi_f,
				       x->fade_in, x->fade_out, x->length,
				       sample_rate, db_to_a(x->amp),
				       x->amp_curve);
  long chan_offset = (long)((out_length * sample_rate) + sweep_length);
  long i;

  for(i=0;i<x->current_num_active_outs;i++)
    x->chan_offset[i] = (i * chan_offset);

  x->T = sweep_length;
  x->T2 = (t_int)(sweep_length + x->chan_offset[x->current_num_active_outs
						     - 1] +
		       (out_length*sample_rate));
  x->current_out_length = out_length;
  x->current_t = 0;
}

void irmeasure_mls_params(t_irmeasure_tilde *x)
{
  double out_length = x->out_length;
  double sample_rate = x->sample_rate;

  long order = x->order;
  long mls_length = (1 << order) - 1;
  long chan_offset = (long) ((out_length*sample_rate) + mls_length);
  long i;

  mls_params(&x->max_length_params, order, db_to_a(x->amp)); /* */

  for(i=0; i<x->current_num_active_outs; i++)
    x->chan_offset[i] = (i * chan_offset);

  x->T = mls_length;
  x->T2 = (t_int) (mls_length + x->chan_offset[x->current_num_active_outs -
						    1] +
			(out_length*sample_rate));
  x->current_out_length = out_length;
  x->current_t = 0;
    
}

void irmeasure_noise_params(t_irmeasure_tilde *x)
{
  t_noise_mode noise_mode = x->noise_mode;

  double length = x->length;
  double out_length = x->out_length;
  double sample_rate = x->sample_rate;

  long noise_length = (long)(sample_rate * length);
  long chan_offset = (long)((out_length * sample_rate) + noise_length);
  long i;

  double amp_comp = 1.;

  if(noise_mode == NOISE_MODE_BROWN)
    amp_comp = x->max_amp_brown;
  if(noise_mode == NOISE_MODE_PINK)
    amp_comp = x->max_amp_pink;

  coloured_noise_params(&x->noise_params, noise_mode, x->fade_in, x->fade_out,
			length, sample_rate, db_to_a(x->amp) / amp_comp); /* */

  for(i=0;i<x->current_num_active_outs;i++)
    x->chan_offset[i] = (i*chan_offset);

  x->T = noise_length;
  x->T2 = (t_int) (noise_length +
			x->chan_offset[x->current_num_active_outs - 1] +
			(out_length * sample_rate));
  x->current_out_length = out_length;
  x->current_t = 0;
}

void irmeasure_params(t_irmeasure_tilde *x)
{
  switch(x->measure_mode)
    {
    case SWEEP:
      irmeasure_sweep_params(x);
      break;
    case MLS:
      irmeasure_mls_params(x);
      break;
    case NOISE:
      irmeasure_noise_params(x);
      break;
    }
}

static double ess_harm_offset(t_ess *x, t_uint harm)
{
  return x->T / log(x->hi_f_act/x->lo_f_act) * log((double) harm);
}

static double max_double(double v1, double v2)
{
  v1 = v1 > v2 ? v1 : v2;
  return v1;
}

static double min_double(double v1, double v2)
{
  v1 = v1 < v2 ? v1 : v2;
  return v1;
}

t_uint32 get_next_lfsr_int(t_uint32 lfsr, t_uint32 feedback_mask)
{
  return (lfsr >> 1) ^ (t_uint32)((0 - (lfsr & 0x1u)) & feedback_mask);
}

void mls_gen_dsp(t_mls *x, t_sample *out, t_uint N)
{
  t_uint32 lfsr = x->lfsr;
  t_uint32 feedback_mask = x->feedback_mask;
  t_uint32 i;

  t_sample amp = (t_sample)x->amp;
  t_sample two_amp = amp * 2.;

  for(i=0; i<N; i++) {
    *out++ = ((lfsr & 0x1u) * two_amp) - amp;
    lfsr = get_next_lfsr_int(lfsr, feedback_mask);
  }

  x->lfsr = lfsr;
}

void mls_gen_block(t_mls *x, t_sample *out, t_uint N)
{
  mls_gen_dsp(x, out, N);
}

void mls_gen(t_mls *x, t_sample *out)
{
  mls_gen_dsp(x, out, x->T);
}

void coloured_noise_gen_dsp(t_noise_params *x, t_sample *out, t_uint startN, t_uint N)
{
  double prev_output = x->prev_output;
  double alpha = x->alpha;

  double alpha0 = x->alpha0;
  double alpha1 = x->alpha1;
  double alpha2 = x->alpha2;
  double alpha3 = x->alpha3;
  double alpha4 = x->alpha4;

  double b0 = x->b0;
  double b1 = x->b1;
  double b2 = x->b2;
  double b3 = x->b3;
  double b4 = x->b4;
  double b5 = x->b5;
  double b6 = x->b6;

  double one_amp = x->amp;
  double two_amp = one_amp * 2.;

  double sample_rate = x->sample_rate;
  double FiN = x->fade_in * sample_rate * 2.;
  double FoN = x->fade_out * sample_rate * 2.;
  double fade_in;
  double fade_out;
  double input_val;
  double result;

  t_uint32 xr = x->gen.x;
  t_uint32 yr = x->gen.y;
  t_uint32 zr = x->gen.z;
  t_uint32 wr = x->gen.w;
  t_uint32 r;

  t_uint T = x->T;
  t_uint i;

  FiN = (FiN < 1.) ? 1. : FiN;
  FoN = (FoN < 1.) ? 1. : FoN;

  switch(x->mode)
    {
    case NOISE_MODE_WHITE:
      for(i=startN;i<startN + N;i++) {
	fade_in = (1-cos(M_PI*min_double(0.5, i/FiN)));
	fade_out = (1-cos(M_PI*min_double(0.5, (T-i)/FoN)));

	r = (xr ^ (xr << 20)) ^ (yr ^ (yr >> 11)) ^ (zr ^ (zr << 27)) ^
	  (wr ^ (wr >> 6));
	xr = yr;
	yr = zr;
	zr = wr;
	wr = r;

	input_val = (r * UNSIGNED_INT32_TO_NORM_DOUBLE * two_amp) - one_amp;
	*out++ = fade_in * fade_out * input_val;
      }
      break;

    case NOISE_MODE_BROWN:
      for(i=startN; i<startN+N; i++) {
	fade_in = (1-cos(M_PI*min_double(0.5, i/FiN)));
	fade_out = (1-cos(M_PI*min_double(0.5, (T-i)/FoN)));
	
	r = (xr ^ (xr << 20)) ^ (yr ^ (yr >> 11)) ^ (zr ^ (zr << 27)) ^
	  (wr ^ (wr >> 6));
	xr = yr;
	yr = zr;
	zr = wr;
	wr = r;

	input_val = ((r * UNSIGNED_INT32_TO_NORM_DOUBLE * two_amp) - one_amp);
	result = prev_output + (alpha * (input_val - prev_output));
       	*out++ = fade_in * fade_out * result;
	prev_output = result;
      }
      break;

    case NOISE_MODE_PINK:
      for(i=startN;i<startN+N;i++) {
	fade_in = (1-cos(M_PI*min_double(0.5,i/FiN)));
	fade_out = (1-cos(M_PI*min_double(0.5,(T-i)/FoN)));
	
	r = (xr ^ (xr << 20)) ^ (yr ^ (yr >> 11)) ^ (zr ^ (zr << 27)) ^
	  (wr ^ (wr >> 6));
	xr = yr;
	yr = zr;
	zr = wr;
	wr = r;

	input_val = (r * UNSIGNED_INT32_TO_NORM_DOUBLE * two_amp) - one_amp;

	b0 = b0 + (alpha0*((input_val * 48.69991228070175) - b0));
	b1 = b1 + (alpha1*((input_val * 11.23890718562874) - b1));
	b2 = b2 + (alpha2*((input_val * 4.96296774193548) - b2));
	b3 = b3 + (alpha3*((input_val * 2.32573483146067) - b3));
	b4 = b4 + (alpha4*((input_val * 1.18433822222222) - b4));
	b5 = -0.7616 * b5 - input_val * 0.0168980;
	result = (b0+b1+b2+b3+b4+b5+b6+input_val*0.5362);
	b6 = input_val * 0.115926;

	*out++ = fade_in * fade_out * result;
	prev_output = result;
      }
      break;
    }

  x->b0 = b0;
  x->b1 = b1;
  x->b2 = b2;
  x->b3 = b3;
  x->b4 = b4;
  x->b5 = b5;
  x->b6 = b6;

  x->gen.x = xr;
  x->gen.y=  yr;
  x->gen.z = zr;
  x->gen.w = wr;

  x->prev_output = prev_output;
}

void coloured_noise_measure(t_noise_params *x, t_uint N, double *max_out_pink, double *max_out_brown)
{
  double prev_output = x->prev_output;
  double alpha = x->alpha;

  double alpha0 = x->alpha0;
  double alpha1 = x->alpha1;
  double alpha2 = x->alpha2;
  double alpha3 = x->alpha3;
  double alpha4 = x->alpha;

  double b0 = x->b0;
  double b1 = x->b1;
  double b2 = x->b2;
  double b3 = x->b3;
  double b4 = x->b4;
  double b5 = x->b5;
  double b6 = x->b6;

  double max_brown = 0.;
  double max_pink = 0.;

  double input_val;
  double result;

  t_uint32 xr = x->gen.x;
  t_uint32 yr = x->gen.y;
  t_uint32 zr = x->gen.z;
  t_uint32 wr = x->gen.w;
  t_uint32 r;

  t_uint i;

  for(i=0;i<N;i++) {
    r = (xr ^ (xr << 20)) ^ (yr ^ (yr >> 11)) ^ (zr ^ (zr << 27)) ^
      (wr ^ (wr >> 6));
    xr = yr;
    yr = zr;
    zr = wr;
    wr = r;

    input_val = ((r * UNSIGNED_INT32_TO_NORM_DOUBLE * 2.) - 1.);
    result = prev_output + (alpha * (input_val - prev_output));
    prev_output = result;

    max_brown = fabs(result) >= max_brown ? fabs(result) : max_brown;
  }

  coloured_noise_reset(x);

  for(i=0;i<N;i++) {
    r = (xr ^ (xr << 20)) ^ (yr ^ (yr >> 11)) ^ (zr ^ (zr << 27)) ^
      (wr ^ (wr >> 6));
    xr = yr;
    yr = zr;
    zr = wr;
    wr = r;

    input_val = (r * UNSIGNED_INT32_TO_NORM_DOUBLE * 2.) - 1.;

    b0 = b0 + (alpha0 * ((input_val * 48.69991228070175) - b0));
    b1 = b1 + (alpha1 * ((input_val * 11.23890718562874) - b1));
    b2 = b2 + (alpha2 * ((input_val * 4.96296774193548) - b2));
    b3 = b3 + (alpha3 * ((input_val * 2.32573483146067) - b3));
    b4 = b4 + (alpha4 * ((input_val * 1.18433822222222) - b4));
    b5 = -0.7616 * b5 - input_val * 0.0168980;
    result = (b0 + b1 + b2 + b3 + b4 + b5 + b6 + input_val * 0.5362);

    b6 = input_val * 0.115926;

    max_pink = fabs(result) >= max_pink ? fabs(result) : max_pink;
  }

  *max_out_pink = max_pink;
  *max_out_brown = max_brown;
}

t_uint ess_igen_dsp(t_ess *x, t_sample *out, t_uint startN, t_uint N, t_bool inv_amp)
{
  double *amp_specifier = x->amp_specifier;

  double K1 = x->K1;
  double K2 = x->K2;
  double amp = x->amp;
  double sample_rate = x->sample_rate;
  double FiN = x->fade_in * sample_rate * 2.;
  double FoN = x->fade_out * sample_rate * 2.;
  double amp_const = (inv_amp == true) ? (4. / amp) * x->lo_f_act * K2 : amp;
  
  double val, fade_in, fade_out, time_val, interp, curve_db, curve_amp;

  t_uint T = x->T;
  t_uint i;
  t_uint j = 2 * x->num_amp_specifiers;

  if(startN>T)
    return 0;

  N = (startN + N > T) ? T - startN : N;
  FiN = (FiN < 1.) ? 1. : FiN;
  FoN = (FoN < 1.) ? 1. : FoN;

  for(i=startN; i<startN+N; i++) {
    fade_in = (1-cos(M_PI*min_double(0.5, (T-i-1)/FiN)));
    fade_out = (1-cos(M_PI*min_double(0.5,(i+1)/FoN)));
    time_val = (T-i-1) * K2;

    for(; time_val<amp_specifier[j]; j-=2);

    interp = (time_val - amp_specifier[j]) / (amp_specifier[j+2] -
					      amp_specifier[j]);
    curve_db = amp_specifier[j+1] + interp*(amp_specifier[j+3] -
					    amp_specifier[j+1]); 
    if(isnan(fabs(curve_db)))
      curve_db = 0.;
    curve_amp = pow(10., -curve_db/20.);

    val = curve_amp * amp_const * fade_in * fade_out * exp(time_val) *
      sin(K1 * (exp(time_val) - 1));
    *out++ = (t_sample)val;
  }

  return N;
}

t_uint ess_igen(t_ess *x, t_sample *out, t_bool inv_amp)
{
  return ess_igen_dsp(x, out, 0, x->T, inv_amp);
}

t_uint ess_gen_dsp(t_ess *x, t_sample *out, t_uint startN, t_uint N)
{
  double *amp_specifier = x->amp_specifier;

  double K1 = x->K1;
  double K2 = x->K2;
  double amp = x->amp;
  double sample_rate = x->sample_rate;
  double FiN = x->fade_in * sample_rate * 2.0;
  double FoN = x->fade_out * sample_rate * 2.0;
  double val, fade_in, fade_out, time_val, interp, curve_db, curve_amp;

  t_uint T = x->T;
  t_uint i;
  t_uint j = 0;

  if(startN > T)
    return 0;

  N = (startN + N > T) ? T - startN : N;
  FiN = (FiN < 1.0) ? 1.0 : FiN;
  FoN = (FoN < 1.0) ? 1.0 : FoN;

  for(i=startN; i<startN+N; i++) {
    fade_in = (1 - cos(M_PI * min_double(0.5, i / FiN)));
    fade_out = (1 - cos(M_PI * min_double(0.5, (T - i) / FoN)));
    time_val = i * K2;

    for(; time_val > amp_specifier[j+2]; j+= 2);
    
    interp = (time_val - amp_specifier[j]) /
      (amp_specifier[j+2] - amp_specifier[j]);
    curve_db = amp_specifier[j+1] +
      interp * (amp_specifier[j+3] - amp_specifier[j+1]);
    curve_amp = pow(10.0, curve_db /  20.);

    val = curve_amp * amp * fade_in * fade_out * sin(K1 * (exp(time_val)-1));
    *out++ = (t_sample)val;
  }
    
  return N;
}

void coloured_noise_gen(t_noise_params *x, t_sample *out)
{
  coloured_noise_gen_dsp(x, out, 0, x->T);
}

void coloured_noise_gen_block(t_noise_params *x, t_sample *out, t_uint startN, t_uint N)
{
  coloured_noise_gen_dsp(x, out, startN, N);
}

t_uint ess_gen(t_ess *x, t_sample *out)
{
  return ess_gen_dsp(x, out, 0, x->T);
}

t_uint ess_gen_block(t_ess *x, t_sample *out, t_uint startN, t_uint N)
{
  return ess_gen_dsp(x, out, startN, N);
}

double irmeasure_param_check(t_irmeasure_tilde *x, char *name, double val,
			  double min, double max)
{
  double new_val = val;
  t_bool changed = false;

  if(val<min) {
    changed = true;
    new_val = min;
  }
  
  if(val>max) {
    changed = true;
    new_val = max;
  }

  if(changed == true)
    pd_error(x,"parameter out of range: setting %s to %lf", name, new_val);

  return new_val;
}

static void irmeasure_tilde_active_ins(t_irmeasure_tilde *x, t_float num_active_ins)
{
  if(num_active_ins < 1) {
    pd_error(x, "at least one input channel must be active");
    num_active_ins = 1;
  }

  if(num_active_ins > x->num_in_chans) {
    pd_error(x, "cannot have more active inputs that actual inputs");
    num_active_ins = x->num_in_chans;
  }

  x->num_active_ins = (t_int)num_active_ins;
}

static void irmeasure_tilde_active_outs(t_irmeasure_tilde *x, t_float num_active_outs)
{
  if(num_active_outs < 1) {
    pd_error(x, "at least one output channel must be active");
    num_active_outs = 1;
  }

  if(num_active_outs > x->num_out_chans) {
    pd_error(x, "cannot have more active outputs that actual outputs");
    num_active_outs = x->num_out_chans;
  }

  x->num_active_outs = (t_int)num_active_outs;
}

static void irmeasure_tilde_extract(t_irmeasure_tilde *x, t_symbol *sym, int argc, t_atom *argv)
{
  t_sample *rec_mem;

  /* int error; */
  long in_chan = 1;
  t_symbol *buffer = NULL;

  t_int rec_length = x->T2;

  t_uint fft_size = x->fft_size;

  rec_mem = x->rec_mem;

  /* Get arguments */

  if (!argc) {
    pd_error(x, "no arguments for extract method");
    return;
  }

  in_chan = argc > 1 ? (long)atom_getint(argv++) : 1;
  buffer = atom_getsymbol(argv++);

  /* Range check */

  if (in_chan < 1 || in_chan > x->current_num_active_ins) {
    pd_error(x, "input channel %ld out of range", in_chan);
    return;
  }

  /* Decrement input channel (reference to zero) */

  in_chan--;

  if (!fft_size) {
    pd_error(x, "no stored impulse responses - you may still be recording");
    return;
  }

  /* Write to buffer */
  buffer_write(buffer, rec_mem + rec_length * in_chan, rec_length, x->resize, 1.0);

}

static void irmeasure_tilde_getir(t_irmeasure_tilde *x, t_symbol *sym, int argc,
			       t_atom *argv)
{
  t_symbol *buffer = NULL;
  
  t_sample *out_buf = NULL;
  t_sample *out_mem = x->out_mem;

  t_uint fft_size = x->fft_size;
  t_uint mem_size = x->out_mem_size;

  t_int T_minus;
  t_int write_length;
  t_int T;

  t_int harmonic = 1;
  t_int in_chan = 1;
  t_int out_chan = 1;

  int err = 0;
  
  if(!argc) {
    pd_error(x, "no arguments for getir message");
    return;
  }

  if(argv->a_type == A_SYMBOL && argc > 1) {
    buffer = atom_getsymbol(argv++);
    in_chan = atom_getint(argv++);
    argc--;
  }
  else 
    buffer = atom_getsymbol(argv++);

  if(argc > 1)
    out_chan = atom_getint(argv++);

  if(argc > 2)
    harmonic = atom_getint(argv++);
  
  if(!fft_size) {
    pd_error(x, "no stored impulse responses - you may still be recording");
    return;
  }

  if(in_chan < 1 || in_chan > x->current_num_active_ins) {
    pd_error(x, "input channel %ld out of range", in_chan); /* */
    return;
  }

  if(out_chan < 1 || out_chan > x->current_num_active_outs) {
    pd_error(x, "output channel %ld out of range", out_chan); /* */
    return;
  }

  if(harmonic != 1 && x->measure_mode != SWEEP) {
    verbose(2, "harmonics only available using sweep measurement");
    harmonic = 1;
  }

  if(harmonic < 1) {
    verbose(2, "harmonic must be an integer of one or greater");
    harmonic = 1;
  }
  
  in_chan--;
  out_chan--;

  if(mem_size < fft_size * x->current_num_active_ins) {
    pd_error(x, "storage memory is not large enough");
    return;
  }

  if(x->measure_mode == SWEEP)
    T_minus = (t_int)ess_harm_offset(&x->sweep_params, harmonic);
  else
    T_minus = 0;

  write_length = (t_int)(x->current_out_length * x->sample_rate);

  if(harmonic > 1) {
    t_int T_minus_prev = (t_int)ess_harm_offset(&x->sweep_params, harmonic - 1);
    t_int L2 = T_minus - T_minus_prev;

    if(L2 < write_length)
      write_length = L2;
  }

  T = x->chan_offset[out_chan] - T_minus;

  while(T < 0)
    T += fft_size;

  out_buf = out_mem + (in_chan * fft_size) + T;
  
  err = buffer_write(buffer, out_buf, write_length, x->resize, 1.);
  
  if(!err)
    pd_error(x, "buffer write error");

  /* buffer_write(buffer,out_buf,write_length,x->write_chan,x->resize,1.); */
  /* rendere vettoriale */
}

void irmeasure_tilde_tone(t_irmeasure_tilde *x, t_float freq, t_float chan)
{
  x->test_tone_freq = freq;

  if(chan < 1)
    chan = -1;

  if(chan > x->num_out_chans)
    chan = x->num_out_chans;

  x->start_measurement = 0;
  x->stop_measurement = 0;
  x->test_tone = (t_int)chan;
}

static void irmeasure_tilde_stop(t_irmeasure_tilde *x)
{
  x->start_measurement = 0;
  x->test_tone = 0;
  x->stop_measurement = 1;
}

static void irmeasure_tilde_clear(t_irmeasure_tilde *x)
{
  irmeasure_tilde_stop(x);

  x->fft_size = 0;
  aligned_freebytes(x->rec_mem, x->rec_mem_size*sizeof(t_sample));
  aligned_freebytes(x->out_mem, x->out_mem_size*sizeof(t_sample));
  x->rec_mem = NULL;
  x->out_mem = NULL;
  x->rec_mem_size = 0;
  x->out_mem_size = 0;
}

static void irmeasure_tilde_noise(t_irmeasure_tilde *x, t_symbol *s, long argc,
				  t_atom *argv)
{
  double length = 10000.;
  double fade_in = 10.;
  double fade_out = 10.;
  double out_length = 5000.;

  t_int num_active_ins = x->num_active_ins;
  t_int num_active_outs = x->num_active_outs;

  t_int mem_size;

  t_noise_mode noise_mode = NOISE_MODE_WHITE;

  if(s == gensym("brown"))
    noise_mode = NOISE_MODE_BROWN;
  if(s == gensym("pink"))
    noise_mode = NOISE_MODE_PINK;

  if(argc > 0)
    length = atom_getfloat(argv++);
  if(argc > 1)
    fade_in = atom_getfloat(argv++);
  if(argc > 2)
    fade_out = atom_getfloat(argv++);
  if(argc > 3)
    out_length = atom_getfloat(argv++);

  length = irmeasure_param_check(x, "length", length, 0., HUGE_VAL);
  fade_in = irmeasure_param_check(x, "fade in time", fade_in, 0., length / 2.);
  fade_out = irmeasure_param_check(x, "fade out time", fade_out, 0.,
				   length / 2.);
  out_length = irmeasure_param_check(x, "ir length", out_length, 0., HUGE_VAL);

  x->noise_mode = noise_mode;
  x->length = length / 1000.;
  x->fade_in = fade_in / 1000.;
  x->fade_out = fade_out / 1000.;
  x->out_length = out_length / 1000.;

  mem_size = num_active_ins * irmeasure_calc_noise_mem_size(x->length,
							    num_active_outs,
							    x->out_length,
							    x->sample_rate);
  if(x->rec_mem_size != mem_size) {
    x->rec_mem = (t_sample *)aligned_resizebytes(x->rec_mem, x->rec_mem_size, mem_size);
    x->rec_mem_size = mem_size;
  }

  x->current_num_active_ins = num_active_ins;
  x->current_num_active_outs = num_active_outs;

  x->measure_mode = NOISE;
  x->fft_size = 0;
  x->test_tone = 0;
  x->stop_measurement = 0;
  x->start_measurement = 1;
}

static void irmeasure_tilde_mls(t_irmeasure_tilde *x, t_symbol *sym, long argc,
			     t_atom *argv)
{
  t_float out_length = 5000.;

  t_int num_active_ins = x->num_active_ins;
  t_int num_active_outs = x->num_active_outs;

  t_int order = 18;

  t_int mem_size;

  if(argc>0)
    order = atom_getint(argv++);
  if(argc>1)
    out_length = atom_getfloat(argv++);

  order = irmeasure_param_check(x, "order", (double)order, 1, 24);
  out_length = irmeasure_param_check(x, "ir length", out_length, 0., HUGE_VAL);

  x->order = (t_int)order;
  x->out_length = out_length / 1000.;

  mem_size = num_active_ins * irmeasure_calc_mls_mem_size(x->order,
							  num_active_outs,
							  x->out_length,
							  x->sample_rate);
  if(x->rec_mem_size != mem_size) {
    x->rec_mem = (t_sample *)aligned_resizebytes(x->rec_mem, x->rec_mem_size, mem_size);
    x->rec_mem_size = mem_size;
  }

  x->current_num_active_ins = num_active_ins;
  x->current_num_active_outs = num_active_outs;

  x->measure_mode = MLS;
  x->fft_size = 0;
  x->test_tone = 0;
  x->stop_measurement = 0;
  x->start_measurement = 1;
}

static void irmeasure_tilde_sweep(t_irmeasure_tilde *x, t_symbol *sym, long argc,
			       t_atom *argv)
{
  t_ess sweep_params;

  double f1 = 20.;
  double f2 = sys_getsr() / 2.;
  double length = 30000.;
  double fade_in = 50.;
  double fade_out = 10.;
  double out_length = 5000.;

  long num_active_ins = x->num_active_ins;
  long num_active_outs = x->num_active_outs;

  t_int mem_size;

  if(argc > 0)
    f1 = atom_getfloat(argv++);
  if(argc > 1)
    f2 = atom_getfloat(argv++);
  if(argc > 2)
    length = atom_getfloat(argv++);
  if(argc > 3)
    fade_in = atom_getfloat(argv++);
  if(argc > 4)
    fade_out = atom_getfloat(argv++);
  if(argc>5)
    out_length = atom_getfloat(argv++);

  f1 = irmeasure_param_check(x, "low frequency", f1, 0.0001,
			     x->sample_rate / 2.);
  f2 = irmeasure_param_check(x, "high frequency", f2, f2, x->sample_rate / 2.);
  length = irmeasure_param_check(x, "length", length, 0., HUGE_VAL);
  fade_in = irmeasure_param_check(x, "fade in time", fade_in, 0., length / 2.);
  fade_out = irmeasure_param_check(x, "fade out time", fade_out, 0.,
				   length / 2.);
  out_length = irmeasure_param_check(x, "ir length", out_length, 0., HUGE_VAL);

  x->lo_f = f1;
  x->hi_f = f2;
  x->fade_in = fade_in / 1000.;
  x->fade_out = fade_out / 1000.;
  x->length = length / 1000.;
  x->out_length = out_length / 1000.;
  
  if(ess_params(&sweep_params, x->lo_f, x->hi_f, x->fade_in, x->fade_out,
		x->length, x->sample_rate, db_to_a(x->amp), NULL)) { /* NULL */
    mem_size = num_active_ins * irmeasure_calc_sweep_mem_size(&sweep_params,
							      num_active_outs,
							      x->out_length,
							      x->sample_rate);
    if(x->rec_mem_size!=mem_size) {
      x->rec_mem = (t_sample *)aligned_resizebytes(x->rec_mem,x->rec_mem_size,mem_size);
      x->rec_mem_size = mem_size;
    }
    if(!x->rec_mem) {
      pd_error(x,"not able to allocate adequate memory for recording");
      x->stop_measurement = 1;
    }
    fill_amp_curve_specifier(x->amp_curve, x->amp_curve_specifier,
			     x->amp_curve_num_specifiers);

    x->current_num_active_ins = num_active_ins;
    x->current_num_active_outs = num_active_outs;

    x->measure_mode = SWEEP;
    x->fft_size = 0;
    x->test_tone = 0;
    x->stop_measurement = 0;
    x->start_measurement = 1;
  }
  else {
    pd_error(x,"zero length sweep");
    x->stop_measurement = 1;
  }
}

static void irmeasure_process(t_irmeasure_tilde *x)
{
  FFT_Setup *fft_setup;

  FFT_Split spectrum_1;
  FFT_Split spectrum_2;
  FFT_Split spectrum_3;
  
  void *measurement_rec;
  void *rec_mem;
  t_sample *excitation_sig = NULL;
  t_sample *out_buf;
  t_sample *out_mem;
  t_float *filter_in;

  t_symbol *filter = filter_retriever(x->deconvolve_filter_specifier);

  t_float filter_specifier[PIRO_MAX_SPECIFIER_ITEMS];
  t_float range_specifier[PIRO_MAX_SPECIFIER_ITEMS];

  t_sample test_pow;
  t_sample max_pow;
  double sample_rate = x->sample_rate;
  double deconvolve_phase = phase_retriever(x->deconvolve_phase);
  /* 0 di default */

  long deconvolve_mode = x->deconvolve_mode;
  long bandlimit = x->measure_mode == SWEEP ? x->bandlimit : 0;

  t_int rec_length = x->T2;
  t_int gen_length = 0;

  t_garray *a;
  t_word *a_filter;

  int filter_length = 0;

  t_uint fft_size;
  t_uint fft_size_log2 = 0;
  t_uint i;

  t_ess sweep_params;
  t_mls max_length_params;
  t_noise_params noise_params;

  t_uint out_mem_size_old;
  
  switch(x->measure_mode)
    {
    case SWEEP:
      ess_params(&sweep_params, x->sweep_params.rf1, x->sweep_params.rf2,
      		 x->sweep_params.fade_in, x->sweep_params.fade_out,
      		 x->sweep_params.RT, x->sweep_params.sample_rate,
      		 x->inv_amp ? x->sweep_params.amp : 1, x->amp_curve);
      gen_length = ess_get_length(&sweep_params);
      break;
    case MLS:
      mls_params(&max_length_params, x->max_length_params.order,
      		 x->inv_amp ? x->max_length_params.amp : 1);
      gen_length = mls_get_length(&max_length_params);
      break;
    case NOISE:
      coloured_noise_params(&noise_params, x->noise_params.mode,
			    x->noise_params.fade_in, x->noise_params.fade_out,
			    x->noise_params.RT, x->noise_params.sample_rate,
			    x->inv_amp ? x->noise_params.amp : 1);
      gen_length = coloured_noise_get_length(&noise_params);
      break;
    }

  fft_size = calculate_fft_size(rec_length + gen_length, &fft_size_log2);
  fft_setup = create_setup(fft_size_log2);

  excitation_sig = (t_sample *)aligned_getbytes(gen_length * sizeof(t_sample));

  spectrum_1.realp = (t_sample *)aligned_getbytes(fft_size * 4 * sizeof(t_sample));
  spectrum_1.imagp = spectrum_1.realp + (fft_size >> 1);
  spectrum_2.realp = spectrum_1.imagp + (fft_size >> 1);
  spectrum_2.imagp = spectrum_2.realp + (fft_size >> 1);
  spectrum_3.realp = spectrum_2.imagp + (fft_size >> 1);
  spectrum_3.imagp = spectrum_3.realp + fft_size;

  attach_array(filter, &a, &a_filter, &filter_length);
  filter_in = filter_length ? (t_float *)aligned_getbytes(fft_size *
							  sizeof(t_float)) : NULL;
  
  if(!fft_setup || !excitation_sig || !spectrum_1.realp || (filter_length &&
							    !filter_in)) {
    pd_error(x, "could not allocate temporary memory");
    destroy_setup(fft_setup);
    aligned_freebytes(excitation_sig, gen_length*sizeof(t_sample));
    aligned_freebytes(spectrum_1.realp, fft_size*4*sizeof(t_sample));
    aligned_freebytes(filter_in,filter_length ? fft_size*sizeof(t_float) : 0);

    return;
  }

  rec_mem = x->rec_mem;
  out_mem_size_old = x->out_mem_size;
  if(out_mem_size_old!=fft_size*x->current_num_active_ins*sizeof(t_sample)) {
    x->out_mem_size = fft_size*x->current_num_active_ins*sizeof(t_sample);
    x->out_mem = (t_sample *)aligned_resizebytes(x->out_mem,out_mem_size_old,
				       x->out_mem_size);
  }
  out_mem = x->out_mem;
  if(!out_mem) {
    pd_error(x,"could not allocate memory");
    aligned_freebytes(excitation_sig, gen_length*sizeof(t_sample));
    aligned_freebytes(filter_in, filter_length ? fft_size*sizeof(t_float) : 0);
    aligned_freebytes(spectrum_1.realp, fft_size*4*sizeof(t_sample));
    destroy_setup(fft_setup);
    
    return;
  }

  switch(x->measure_mode)
    {
    case SWEEP:
      ess_gen(&sweep_params, excitation_sig);
      break;
    case MLS:
      mls_gen(&max_length_params, excitation_sig);
      break;
    case NOISE:
      coloured_noise_gen(&noise_params, excitation_sig);
      break;
    }
  
  time_to_halfspectrum(fft_setup, excitation_sig, gen_length, spectrum_2,
		       fft_size);
  
  if(bandlimit) {
    ess_igen(&sweep_params, excitation_sig, true);
    
    time_to_halfspectrum(fft_setup, excitation_sig, gen_length, spectrum_3,
			 fft_size);
    convolve(spectrum_3, spectrum_2, fft_size, SPECTRUM_REAL);
    power_full_spectrum_from_half_spectrum(spectrum_3, fft_size);
    variable_phase_from_power_spectrum(fft_setup, spectrum_3, fft_size,
    				       deconvolve_phase, true);

    spectrum_3.imagp[0] = spectrum_3.realp[fft_size >> 1];
  }
  else {
    for(i=1, max_pow=0.; i<(fft_size>>1); i++) {
      test_pow = spectrum_2.realp[i]*spectrum_2.realp[i] +
      	spectrum_2.imagp[i]*spectrum_2.imagp[i];
      max_pow = test_pow > max_pow ? test_pow : max_pow;
    }
    
    max_pow = pow_to_db(max_pow);
    
    fill_power_array_specifier(filter_specifier, x->deconvolve_filter_specifier,
  			       x->deconvolve_num_filter_specifiers);    
    fill_power_array_specifier(range_specifier, x->deconvolve_range_specifier,
  			       x->deconvolve_num_range_specifiers);    
    buffer_read(filter, filter_in, fft_size);
    make_deconvolution_filter(fft_setup, spectrum_2, spectrum_3,
    			      filter_specifier, range_specifier, max_pow,
    			      filter_in, fft_size, fft_size,
    			      SPECTRUM_REAL, deconvolve_mode, deconvolve_phase,
    			      sample_rate);   
  }

  for(i=0; i<(t_uint)x->current_num_active_ins; i++) {
    measurement_rec = (t_sample *)rec_mem + (i * rec_length);
    out_buf = out_mem + (i * fft_size); 
    
    time_to_halfspectrum(fft_setup, measurement_rec, rec_length, spectrum_1,
			 fft_size);
    
    deconvolve_with_filter(spectrum_1, spectrum_2, spectrum_3, fft_size,
    			   SPECTRUM_REAL);
    
    spectrum_to_time(fft_setup, out_buf, spectrum_1, fft_size, SPECTRUM_REAL); 
    
  }

  destroy_setup(fft_setup);
  aligned_freebytes(excitation_sig, gen_length*sizeof(t_sample));
  aligned_freebytes(spectrum_1.realp, fft_size*4*sizeof(t_sample));
  aligned_freebytes(filter_in, filter_length ? fft_size*sizeof(t_float) : 0);
  
  x->fft_size = fft_size;
  
  outlet_bang(x->bangout);
}

static void irmeasure_tilde_reprocess(t_irmeasure_tilde *x)
{
  if(!x->fft_size) {
    pd_error(x, "no complete recording to process");
    return;
  }
  
  irmeasure_process(x);
}

static void irmeasure_tilde_inv_amp(t_irmeasure_tilde *x, t_float f)
{
  t_uint i = (t_uint)f;

  x->inv_amp = f ? 1 : 0;
}

static void irmeasure_tilde_bandlimit(t_irmeasure_tilde *x, t_float f)
{
  t_uint i = (t_uint)f;

  x->bandlimit = (f >= 0 && f <=1 ) ? f : 1;
}

static void irmeasure_tilde_fullscale(t_irmeasure_tilde *x, t_float f)
{
  x->amp = f; /* non metto check perche' inclusi in db_to_a */
}

static void irmeasure_tilde_deconvdelay(t_irmeasure_tilde *x, t_symbol *sym,
					long argc, t_atom *argv)
{
  if(argc && argv && (argv->a_type == A_FLOAT)) {
    x->deconvolve_delay = *argv;
  }
  else
    SETSYMBOL(&x->deconvolve_delay, gensym("center"));

  if(argc && argv && argv->a_type == A_SYMBOL && atom_getsymbol(argv) != gensym("center"))
    pd_error(x, "unknown symbol for the deconvolution delay (using center)");
}

static void irmeasure_tilde_deconvphase(t_irmeasure_tilde *x, t_symbol *sym,
					long argc, t_atom *argv)
{
  if(argc && argv)
    x->deconvolve_phase = phase_parser(*argv);
  else
    SETSYMBOL(&x->deconvolve_phase, gensym("minimum"));

}

static void irmeasure_tilde_deconvfilter(t_irmeasure_tilde *x, t_symbol *sym,
				      long argc, t_atom *argv)
{
  long i;

  if(!argv) {
    SETFLOAT(x->deconvolve_filter_specifier, 20.);
    SETFLOAT(x->deconvolve_filter_specifier+1, -20.);
    SETFLOAT(x->deconvolve_filter_specifier+2, 30.);
    SETFLOAT(x->deconvolve_filter_specifier+3, -100.);
    SETFLOAT(x->deconvolve_filter_specifier+4, 19000.);
    SETFLOAT(x->deconvolve_filter_specifier+5, -100.);
    SETFLOAT(x->deconvolve_filter_specifier+6, 22050.);
    SETFLOAT(x->deconvolve_filter_specifier+7, -20.);

    x->deconvolve_num_filter_specifiers = 8;

    return;
  }

  if(argc && argv) {
    if(argc > 1 && argc & 1) {
      pd_error(x, "frequency without level found in filter specification list");
      argc--;
    }

    for(i=0; i<argc; i++, argv++) {
      if(argv->a_type == A_FLOAT || !i) 
	if(!i && argv->a_type == A_SYMBOL) {
	  SETSYMBOL(x->deconvolve_filter_specifier, argv->a_w.w_symbol);
	}
	else 
	  x->deconvolve_filter_specifier[i] = *argv;	
      else {
	SETFLOAT(x->deconvolve_filter_specifier+i, 0);
	pd_error(x, "symbol found in filter specification list");
      }
    }

    x->deconvolve_num_filter_specifiers = argc;
  }
  else {
    SETFLOAT(x->deconvolve_filter_specifier, -100);
    x->deconvolve_num_filter_specifiers = 1;
  }  
}

static void irmeasure_tilde_deconvrange(t_irmeasure_tilde *x, t_symbol *sym,
					long argc, t_atom *argv)
{
  long i;

  if(argc && argv) {
    if(argc > 1 && argc & 1) {
      pd_error(x, "frequency without level found in range specification list");
      argc--;
    }

    for(i=0; i<argc; i++, argv++) {
      if(argv->a_type == A_FLOAT || !i)
	x->deconvolve_range_specifier[i] = *argv;
      else {
	SETFLOAT(x->deconvolve_range_specifier+i, 0);
	pd_error(x, "symbol found in range specification list");
      }
    }

    x->deconvolve_num_range_specifiers = argc;
  }
  else {
    SETFLOAT(x->deconvolve_range_specifier, 1000);
    x->deconvolve_num_range_specifiers = 1;
  }
  
}

static void irmeasure_tilde_deconvmode(t_irmeasure_tilde *x, t_float f)
{
  t_int i = (t_uint)f;

  x->deconvolve_mode = (i >= 0 && i <= 2 ) ? i : 0;
}

static void irmeasure_perform_excitation(t_irmeasure_tilde *x, void *out,
				      long current_t, long vec_size)
{
  long start;
  long todo;

  start = (long)max_double(0, -current_t);

  out = ((t_sample *) out) + start;    

  todo = (long)max_double(0,min_double(vec_size - start,
				       min_double((double)x->T,
						  (double)(x->T - current_t))));

  if(current_t <= 0 && todo) {
    mls_reset(&x->max_length_params); /* */
    coloured_noise_reset(&x->noise_params);
  }

  if(start<vec_size && todo) {
    switch(x->measure_mode)
      {
      case SWEEP:
	ess_gen_block(&x->sweep_params, out, current_t, todo);
	break;
      case MLS:
	mls_gen_block(&x->max_length_params, out, todo);
	break;
      case NOISE:
	coloured_noise_gen_block(&x->noise_params, out, current_t, todo);
	break;
      }
  }
}

static t_int *irmeasure_tilde_perform(t_int *w)
{
  t_irmeasure_tilde *x = (t_irmeasure_tilde *)(w[1]);
  long vec_size = (long)(w[2]);

  t_sample *measurement_rec;
  t_sample *rec_mem;
  t_sample *in;
  t_sample *out;

  double phase = x->phase;
  double phase_inc = x->test_tone_freq / x->sample_rate;
  double amp = db_to_a(x->amp);
  double sample_rate = x->sample_rate;
  double progress_mul = 0.;

  t_int T2;
  t_int current_t;
  t_int current_t2;
  t_int mem_size = 0;
  t_int i, j;

  long test_tone = x->test_tone;
  long excitation_playing;
  long mem_check;
  long current_num_active_ins = x->current_num_active_ins;
  long current_num_active_outs = x->current_num_active_outs;
  long num_out_chans = x->num_out_chans;
  
  if(x->start_measurement)
    irmeasure_params(x);

  if(x->stop_measurement)
    x->current_t = x->T2;

  x->start_measurement = 0;
  x->stop_measurement = 0;
  T2 = x->T2;

  current_t = x->current_t;
  current_t2 = current_t;
  excitation_playing = current_t2 < T2;

  rec_mem = x->rec_mem;

  mem_size = irmeasure_calc_mem_size(x, current_num_active_ins,
				     current_num_active_outs, sample_rate);
  mem_check = x->rec_mem_size >= (t_int)mem_size;

  if(mem_check) {

    /* input */
    
    for(j=0; j<current_num_active_ins; j++) {
      in = (t_sample *)x->in_chans[j];
      measurement_rec = rec_mem + (j * T2);
      current_t2 = current_t;
      for(i=0; i<vec_size && current_t2 < T2; i++, current_t2++)
	measurement_rec[current_t2] = *in++;
    }
  }

  for(j=0; j<num_out_chans+1; j++)
    for(i=0, out=x->out_chans[j]; i<vec_size; i++)
      *out++ = 0.;

  if(test_tone) {
    if(test_tone == -1)
      out = x->out_chans[0];
    else
      out = x->out_chans[test_tone - 1];

    for(i=0; i<vec_size; i++) {
      *out++ = (t_sample )(amp * sin(phase * M_PI * 2));
      phase += phase_inc;
    }

    while(phase<0.)
      phase += 1.;
    while(phase>1.)
      phase -= 1.;

    if(test_tone == -1) {
      for(j=1, in=x->out_chans[0]; j<current_num_active_outs; j++) {
	for(i=0, out=x->out_chans[j]; i<vec_size; i++)
	  out[i] = in[i];
      }
    }
  }
  else {
    if(mem_check) {

      /* output */
      
      for(j=0; j<current_num_active_outs; j++)
	irmeasure_perform_excitation(x, x->out_chans[j],
				     (long)(x->current_t - x->chan_offset[j]),
				     vec_size);
      if(x->abs_progress)
	progress_mul = 1000. / sample_rate;
      else
	if(T2)
	  progress_mul = 1. / T2;

      out = x->out_chans[num_out_chans];
      current_t2 = x->current_t;
      for(i=0; i<vec_size && current_t2 < T2; i++, current_t2++)
	out[i] = (t_sample)(progress_mul * current_t2);
      for(; i<vec_size; i++)
	  out[i] = (t_sample)(progress_mul * T2);
    }
  }

  x->current_t = current_t2;
  x->phase = phase;

  if(excitation_playing && current_t2 >= T2) 
    irmeasure_process(x); /* in Max usa defer: portare con clock_delay(0)... */
  
  return (w+3);
}

static void irmeasure_tilde_dsp(t_irmeasure_tilde *x, t_signal **sp)
{
  long i;
  t_ess sweep_params;
  t_noise_params noise_params;
  double old_sr;

  t_int mem_size = 0;
  t_int old_mem_size = x->rec_mem_size;

  /* store sample rate */

  old_sr = x->sample_rate;
  x->sample_rate = sys_getsr();

  if(x->sample_rate != old_sr || x->no_dsp) {
    coloured_noise_params(&noise_params, 0, 0, 0, 1, x->sample_rate, 1);
    coloured_noise_measure(&noise_params, (1 << 25), &x->max_amp_pink,
			   &x->max_amp_brown);
    x->no_dsp = 0;
  }
  
  if(x->start_measurement) {
    switch(x->measure_mode) {
    case SWEEP:
      ess_params(&sweep_params, x->lo_f, x->hi_f, x->fade_in, x->fade_out,
		 x->length, x->sample_rate, db_to_a(x->amp), 0);
      mem_size = x->current_num_active_ins *
	irmeasure_calc_sweep_mem_size(&sweep_params, x->current_num_active_outs,
				   x->out_length, x->sample_rate);
      break;
    case MLS:
      mem_size = x->current_num_active_ins *
	irmeasure_calc_mls_mem_size(x->order, x->current_num_active_outs,
				 x->out_length, x->sample_rate);
      break;
    case NOISE:
      mem_size = x->current_num_active_ins *
	irmeasure_calc_noise_mem_size(x->length, x->current_num_active_outs,
				   x->out_length, x->sample_rate);
      break;
    }
    if(old_mem_size!=mem_size)
      x->rec_mem = (t_sample *)aligned_resizebytes(x->rec_mem, old_mem_size, mem_size);
    if(!x->rec_mem) {
      mem_size = old_mem_size;
      pd_error(x, "not able to allocate adequate memory for recording");
    }
  }
  else 
    if(x->sample_rate != old_sr)
      irmeasure_tilde_clear(x);

  for(i=0;i<x->num_in_chans;i++)
    x->in_chans[i] = sp[i]->s_vec;
  for(;i<PIRO_MAX_MEASURE_CHANS;i++)
    x->in_chans[i] = 0;

  for(i=0;i<x->num_out_chans;i++)
    x->out_chans[i] = sp[i+x->num_in_chans]->s_vec;
  x->out_chans[i] = sp[i+x->num_in_chans]->s_vec;
  for(i++;i<PIRO_MAX_MEASURE_CHANS+1;i++)
    x->out_chans[i] = 0;
  
  dsp_add(irmeasure_tilde_perform, 2, x, sp[0]->s_n);
}

static void irmeasure_tilde_free(t_irmeasure_tilde *x)
{
  aligned_freebytes(x->rec_mem, x->rec_mem_size);
  aligned_freebytes(x->out_mem, x->out_mem_size);
  aligned_freebytes(x->deconvolve_filter_specifier,
	    PIRO_MAX_SPECIFIER_ITEMS*sizeof(t_atom));
  aligned_freebytes(x->deconvolve_range_specifier,
	    PIRO_MAX_SPECIFIER_ITEMS*sizeof(t_atom));

  outlet_free(x->bangout);
}

static void *irmeasure_tilde_new(t_symbol *s, short argc, t_atom *argv)
{
  t_irmeasure_tilde *x = (t_irmeasure_tilde *)pd_new(irmeasure_tilde_class);

  t_int num_out_chans = 1;
  t_int num_in_chans = 1;
  long i = 0;

  if(argc && argv->a_type == A_FLOAT) {
    num_in_chans = atom_getint(argv++);
    num_in_chans = num_in_chans < 1 ? 1 : num_in_chans;
    num_in_chans = num_in_chans > PIRO_MAX_MEASURE_CHANS ?
      PIRO_MAX_MEASURE_CHANS : num_in_chans;
    argc--;
  }

  if(argc && argv->a_type == A_FLOAT) {
    num_out_chans = atom_getint(argv++);
    num_out_chans = num_out_chans < 1 ? 1 : num_out_chans;
    num_out_chans = num_out_chans > PIRO_MAX_MEASURE_CHANS ?
      PIRO_MAX_MEASURE_CHANS : num_out_chans;
  }

  /* attributes */
  x->bandlimit = 1; /* bandlimit 1 di default */
  x->abs_progress = 0;
  x->amp = 0.; /* -1. in HIRT e aggiunto metodo fullscale */
  x->inv_amp = 0;

  /* deconvolution attributes & methods */
  x->deconvolve_mode = 0; /* FILTER_REGULARISATION */
  x->deconvolve_filter_specifier = (t_atom *)aligned_getbytes(PIRO_MAX_SPECIFIER_ITEMS*sizeof(t_atom));
  x->deconvolve_range_specifier = (t_atom *)aligned_getbytes(PIRO_MAX_SPECIFIER_ITEMS*sizeof(t_atom));
  if(!x->deconvolve_filter_specifier || !x->deconvolve_range_specifier)
    pd_error(x,"could not allocate space for attribute storage");
  irmeasure_tilde_deconvrange(x, 0, 0, 0);
  irmeasure_tilde_deconvfilter(x, 0, 0, 0);
  irmeasure_tilde_deconvphase(x, 0, 0, 0);
  irmeasure_tilde_deconvdelay(x, 0, 0, 0); // sicuro?

  x->amp_curve_num_specifiers = 0;

  x->write_chan = 1;
  x->resize = 1; /* faccio sempre il resize */
  
  x->T2 = 0; /* */
  x->current_t = 0;
  x->fft_size = 0;
  x->start_measurement = 0;
  x->stop_measurement = 0;
  x->test_tone = 0;
  x->no_dsp = 1; /* fissata la compensazione in ampiezza dei rumori */
  x->sample_rate = sys_getsr();

  if(!x->sample_rate)
    x->sample_rate = 44100.0;

  x->rec_mem_size = 0;
  x->out_mem_size = 0;
  x->rec_mem = (t_sample *)aligned_getbytes(x->rec_mem_size*sizeof(t_sample));
  x->out_mem = (t_sample *)aligned_getbytes(x->out_mem_size*sizeof(t_sample));
  if(!x->rec_mem || !x->out_mem)
    pd_error(x,"could not allocate space for recording storage");

  x->num_in_chans = (long)num_in_chans;
  x->num_out_chans = (long)num_out_chans;
  x->num_active_ins = (long)num_in_chans;
  x->num_active_outs = (long)num_out_chans;
  x->current_num_active_ins = (long)num_in_chans;
  x->current_num_active_outs = (long)num_out_chans;

  x->measure_mode = SWEEP; /* default configuration */
  x->phase = 0.;

  /* questa memoria va liberata */
  for(i=1;i<x->num_in_chans;i++)
    inlet_new(&x->x_obj,&x->x_obj.ob_pd,&s_signal,&s_signal);
  
  for(i=0;i<x->num_out_chans+1;i++) /* */
    outlet_new(&x->x_obj, &s_signal);
  /* un outlet e' per abs_progress */

  x->bangout = outlet_new(&x->x_obj, &s_bang);

  return (void *)x;
}

void irmeasure_tilde_setup(void)
{
#if PD_FLOATSIZE == 32
  irmeasure_tilde_class = class_new(gensym("irmeasure~"),
				 (t_newmethod)irmeasure_tilde_new,
				 (t_method)irmeasure_tilde_free,
				 sizeof(t_irmeasure_tilde),
				 CLASS_DEFAULT, 
				 A_GIMME, 0);
#elif PD_FLOATSIZE == 64
  irmeasure_tilde_class = class_new64(gensym("irmeasure~"),
				   (t_newmethod)irmeasure_tilde_new,
				   (t_method)irmeasure_tilde_free,
				   sizeof(t_irmeasure_tilde),
				   CLASS_DEFAULT, 
				   A_GIMME, 0);
#else
#error [irmeasure~]: invalid FLOATSIZE: must be 32 or 64
#endif
  
  CLASS_MAINSIGNALIN(irmeasure_tilde_class, t_irmeasure_tilde, x_f);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_dsp, gensym("dsp"), 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_sweep, gensym("sweep"), A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_mls, gensym("mls"), A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_noise, gensym("white"), A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_noise, gensym("brown"), A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_noise, gensym("pink"), A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_reprocess, gensym("reprocess"), 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_stop, gensym("stop"), 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_clear, gensym("clear"), 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_tone, gensym("tone"), A_FLOAT,
		  A_FLOAT, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_active_ins, gensym("activeins"),
		  A_FLOAT, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_active_outs, gensym("activeouts"),
		  A_FLOAT, 0);
  
  class_addmethod(irmeasure_tilde_class, (t_method)irmeasure_tilde_getir,
		  gensym("getir"), A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class, (t_method)irmeasure_tilde_extract,
		  gensym("extract"), A_GIMME, 0);
  
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_inv_amp, gensym("invamp"),
		  A_FLOAT, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_bandlimit, gensym("bandlimit"),
		  A_FLOAT, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_fullscale, gensym("fullscale"),
		  A_FLOAT, 0);
  
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_deconvmode, gensym("deconvmode"),
		  A_FLOAT, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_deconvrange, gensym("deconvrange"),
		  A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_deconvfilter,
		  gensym("deconvfilter"),A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_deconvphase, gensym("deconvphase"),
		  A_GIMME, 0);
  class_addmethod(irmeasure_tilde_class,
		  (t_method)irmeasure_tilde_deconvdelay, gensym("deconvdelay"),
		  A_GIMME, 0);
}
