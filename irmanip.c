#include "irmanip.h"



static t_class *irmanip_class;

typedef struct _irmanip_tilde {
  t_object  x_obj;
  //t_sample x_f;

  double sample_rate;
	
  /* deconvolution attr */
  t_atom *deconvolve_filter_specifier;
  t_atom *deconvolve_range_specifier;
  long deconvolve_num_filter_specifiers;
  long deconvolve_num_range_specifiers;
  long deconvolve_mode;
  t_atom deconvolve_phase; //
  t_atom deconvolve_delay;

  /* buffer write attr */
  long resize;
  long write_chan;

  t_outlet *floatout;
  
} t_irmanip;

typedef enum {
  RMS_RESULT_SUCCESS = 0,
  RMS_RESULT_IN_LEVEL_NOT_FOUND = 1,
  RMS_RESULT_OUT_LEVEL_NOT_FOUND = 2
} t_rms_result;

typedef struct _rms_measure
{
  t_float accum;
  t_uint last_index;
} t_rms_measure;

void fill_amp_curve_specifier(double *array, t_atom *specifier, long num_specifiers)
{
  long i=0;

  for(i=0;i<num_specifiers;i++) {
    array[i] = atom_getfloat(specifier+i);
  }

  array[i] = HUGE_VAL;
}

t_float norm_find_max(t_float *in, t_uint length, t_float start_max)
{
  t_float norm_factor = start_max;
  t_float norm_test;

  t_uint i;

  for(i=0; i<length; i++) {
    norm_test = fabs(in[i]);

    if(norm_test > norm_factor)
      norm_factor = norm_test;
  }

  return norm_factor;
}

t_float calculate_rms_run(t_rms_measure *rms, t_sample *in, t_uint length, t_uint width, t_uint index)
{
  t_float accum = rms->accum;
  t_float in_val;

  t_uint last_index = rms->last_index;
  t_uint half_length = width >> 1;
  t_uint pre_length = half_length;
  t_uint post_length = half_length + 1;
  t_uint i;

  width = (half_length << 1) + 1;

  if(last_index == index+1) {
    if(index>=half_length) {
      in_val = in[index - half_length];
      accum += in_val * in_val;
    }

    if(index + half_length + 1 < length) {
      in_val = in[index+half_length+1];
      accum -= in_val * in_val;
    }

    rms->accum = accum;
    rms->last_index = index;

    return sqrt(accum / width);
    
  }

  if(last_index == index - 1) {
    if(index >= half_length + 1) {
      in_val = in[index - (half_length + 1)];
      accum -= in_val * in_val;
    }

    if(index + half_length < length) {
      in_val = in[index + half_length];
      accum += in_val * in_val;
    }

    rms->accum = accum;
    rms->last_index = index;

    return sqrt(accum / width);
  }

  pre_length = pre_length > index ? index : pre_length;
  post_length = (index + post_length >= length) ? (length - 1) - index: post_length;

  for(i=index-pre_length, accum=0.; i<=post_length; i++) {
    in_val = in[i];
    accum += in_val * in_val;
  }

  rms->accum = accum;
  rms->last_index = index;

  return sqrt(accum / width);
  
}

void reset_rms(t_rms_measure *rms, t_uint index_reset)
{
  rms->accum = 0.;
  rms->last_index = index_reset;
}

t_rms_result trim_find_crossings_rms(t_sample *in_buf, t_uint length, t_uint window_in, t_uint window_out, t_float in_db, t_float out_db, t_float mul, t_uint *current_start, t_uint *current_end)
{
  t_uint start_search = *current_start;
  t_uint end_search = *current_end;
  t_uint i, j;

  t_float in_lin = pow(10., in_db / 20.);
  t_float out_lin = pow(10., out_db / 20.);

  t_rms_measure rms;

  in_lin = (in_db == -HUGE_VAL) ? -1. : in_lin;
  out_lin = (out_db == -HUGE_VAL) ? -1. : out_lin;

  reset_rms(&rms, length + 2);

  for(i=0; i < length && i < start_search; i++) 
    if(mul * calculate_rms_run(&rms, in_buf, length, window_in, i) >in_lin)
      break;

  if(i==length)
    return RMS_RESULT_IN_LEVEL_NOT_FOUND;

  reset_rms(&rms, length + 2);

  for(j=length; j>i && j>end_search; j--)
    if(mul * calculate_rms_run(&rms, in_buf, length, window_out, j-1) > out_lin)
      break;

  if(j==i)
    return RMS_RESULT_OUT_LEVEL_NOT_FOUND;

  *current_start = i;
  *current_end = j + 1;
      
    
  return RMS_RESULT_SUCCESS;
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

static void irmanip_deconvdelay(t_irmanip *x, t_symbol *sym,
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

static void irmanip_deconvphase(t_irmanip *x, t_symbol *sym,
					long argc, t_atom *argv)
{
  if(argc && argv)
    x->deconvolve_phase = phase_parser(*argv);
  else
    SETSYMBOL(&x->deconvolve_phase, gensym("minimum"));

}

static void irmanip_deconvfilter(t_irmanip *x, t_symbol *sym,
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

static void irmanip_deconvrange(t_irmanip *x, t_symbol *sym,
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

static void irmanip_deconvmode(t_irmanip *x, t_float f)
{
  t_int i = (t_uint)f;

  x->deconvolve_mode = (i >= 0 && i <= 2 ) ? i : 0;
}

int trim_check_number(t_atom *a)
{
  if(a->a_type == A_FLOAT)
    return 0;
  else
    return 1;
}

int trim_check_db(t_atom *a, t_float *db)
{
  if(!trim_check_number(a)) {
    *db = atom_getfloat(a);
    return 0;
  }

  if(a->a_type == A_SYMBOL) {
    if(atom_getsymbol(a) == gensym("off")) {
      *db = -HUGE_VAL;
      return 0;
    }
  }

  return 1;
}

t_float calculate_norm(t_float **samples, t_int *lengths, short N, t_float normalize, int norm_mode)
{
  t_float max = 0.;
  short i;

  if(norm_mode)
    for(i=0; i<N; i++)
      max = norm_find_max(samples[i], lengths[i], max);

  if(max)
    return db_to_a(normalize) / max;
  else
    return 1.;
}

t_float calculate_trim(t_float **samples, t_int *lengths, t_int max_length, t_int total_length, short N, t_float in_db, t_float out_db, t_float sample_rate, t_int *trim_offset, t_int *trim_length, t_float normalize, int norm_mode)
{
  t_uint current_start = max_length;
  t_uint current_end = 0;
  short no_success = 1;
  short i;

  t_float integration_times[2];
  t_float mul;

  integration_times[0] = 5.;
  integration_times[1] = 50.;
  
  mul = calculate_norm(samples, lengths, N, normalize, norm_mode);

  for(i=0; i<N; i++) {
    if(!trim_find_crossings_rms(samples[i], lengths[i], (t_uint)(integration_times[0] * sample_rate / 1000.), (t_uint)(integration_times[1] * sample_rate / 1000.), in_db, out_db, mul, &current_start, &current_end))
      no_success = 0;
  }

  if(no_success) {
    *trim_offset = 0;
    *trim_length = 0;
  }
  else {
    *trim_offset = current_start;
    *trim_length = current_end - current_start;
  }

  return mul;
}

typedef enum {
  FADE_LIN = 0,
  FADE_SQUARE = 1,
  FADE_SQUARE_ROOT = 2,
  FADE_COS = 3,
  FADE_GOMPERTZ = 4,
} t_fade_type;

void fade_calc_fade_out(t_float *in_buf, t_uint fade_length, t_uint length, t_fade_type fade_type)
{
  t_float mult = 1. / fade_length;
  t_float fade_val;

  t_uint i;
  
  if(fade_length > length)
    fade_length = length;

  in_buf += length - 1;

  switch(fade_type)
    {
    case FADE_LIN:
      for(i=0; i<fade_length; i++) {
	fade_val = i * mult;
	*in_buf-- *= fade_val;
      }
      break;

    case FADE_SQUARE:
      for(i=0; i<fade_length; i++) {
	fade_val = i * mult;
	fade_val *= fade_val;
	*in_buf-- *= fade_val;
      }
      break;

    case FADE_SQUARE_ROOT:
      for(i=0; i<fade_length; i++) {
	fade_val = i * mult;
	fade_val = sqrt(fade_val);
	*in_buf-- *= fade_val;
      }
      break;

    case FADE_COS:
      mult *= M_PI * 0.5;
      for(i=0; i< fade_length; i++) {
	fade_val = i * mult;
	fade_val = sin(fade_val);
	*in_buf--  *= fade_val;
      }
      break;

    case FADE_GOMPERTZ:
      mult *= 2.9;
      for(i=0; i<fade_length; i++) {
	fade_val = (i * mult) - 1.2;
	fade_val = exp(-0.1 * exp(-5. * fade_val));
	*in_buf-- *= fade_val;
      }
      break;
      
    }
}

void fade_calc_fade_in(t_float *in_buf, t_uint fade_length, t_uint length, t_fade_type fade_type)
{
  t_float mult = 1. / fade_length;
  t_float fade_val;

  t_uint i;

  if(fade_length > length)
    fade_length = length;

  switch(fade_type)
    {
    case FADE_LIN:
      for(i=0; i<fade_length; i++) {
	fade_val = i * mult;
	*in_buf++ *= fade_val;
      }
      break;

    case FADE_SQUARE:
      for(i=0; i<fade_length; i++) {
	fade_val = i * mult;
	fade_val *= fade_val;
	*in_buf++ *= fade_val;
      }
      break;

    case FADE_SQUARE_ROOT:
      for(i=0; i<fade_length; i++) {
	fade_val = i *mult;
	fade_val = sqrt(fade_val);
	*in_buf++ *= fade_val;
      }
      break;

    case FADE_COS:
      mult *= M_PI * 0.5;
      for(i=0; i<fade_length; i++) {
	fade_val = i * mult;
	fade_val = sin(fade_val);
	*in_buf++ *= fade_val;
      }
      break;

    case FADE_GOMPERTZ:
      mult *= 2.9;
      for(i=0; i<fade_length; i++) {
	fade_val = (i * mult) - 1.2;
	fade_val = exp(-0.1 * exp(-5. * fade_val));
	*in_buf++ *= fade_val;
      }
      break;
      
    }
}

void trim_copy_part(t_float *out_buf, t_float *in_buf, t_uint offset, t_uint length)
{
  t_uint i;

  in_buf += offset;

  for(i=0; i<length; i++)
    *out_buf++ = *in_buf++;
}

void trim_write_internal_buffer(t_float *samples, t_float *internal_buffer, t_int offset, t_int length, t_int fade_in, t_int fade_out, t_int pad_in, t_int pad_out, t_fade_type fade_type)
{
  t_int i;

  for(i=0; i<pad_in; i++)
    internal_buffer[i] = 0.;

  trim_copy_part(internal_buffer+pad_in, samples, offset, length);
  fade_calc_fade_in(internal_buffer+pad_in, fade_in, length, fade_type);
  fade_calc_fade_out(internal_buffer+pad_in, fade_out, length, fade_type);

  for(i=pad_in+length; i<(pad_in+length+pad_out); i++)
    internal_buffer[i] = 0.;
}

long trim_write_buffer(t_symbol *buffer, t_float *samples, t_float *temp_buf, t_int trim_offset, t_int trim_length, t_int fade_in, t_int fade_out, t_int pad_in, t_int pad_out, t_float norm_factor, t_int L, t_float sample_rate, t_fade_type fade_type, t_int end)
{
  long err = 0;
  long length = trim_length;
  if(L < trim_offset + trim_length)
    trim_length = L - trim_offset;

  if(L < trim_offset)
    return err;

  if(trim_offset < fade_in) 
    fade_in = trim_offset;

  if(L < trim_offset + trim_length + fade_out)
    fade_out = L - (trim_offset + trim_length);

  trim_offset -= fade_in;
  trim_length += fade_in + fade_out;

  if(end < trim_length)
    length = end;

  trim_write_internal_buffer(samples, temp_buf, trim_offset, length, fade_in, fade_out, pad_in, pad_out, fade_type);
  
  err = buffer_write(buffer, temp_buf, length+pad_in+pad_out, 1, norm_factor);

  return err;
}

static void irmanip_trim(t_irmanip *x, t_symbol *sym,long argc,
			      t_atom *argv)
{
  t_symbol *in_buffer_names[128];
  t_symbol *out_buffer_names[128];

  t_float *samples[128];
  t_int lengths[128];

  t_float in_db = -HUGE_VAL;
  t_float out_db = -HUGE_VAL;
  t_float fade_in_time;
  t_float fade_out_time;
  t_float pad_in_time = 0.;
  t_float pad_out_time = 0.;
  t_float sample_rate = x->sample_rate;
  t_float norm_factor;

  t_float normalize = 0.;
  int fade_type = 0;
  int norm_mode = 0;

  int argc_halved;

  t_float *temp_buf;

  t_int trim_offset;
  t_int trim_length;
  t_int fade_in;
  t_int fade_out;
  t_int pad_in;
  t_int pad_out;
  t_int overall_length = 0;
  t_int max_length = 0;
  t_int offset = 0;
  t_int i, j;

  t_int end = (t_int)HUGE_VAL;

  short num_buffers = 0;

  t_int in_place = 0; // trim or trimto

  /* t_bool overall_error = false;  */

  if(argc > 4 && !trim_check_number(argv + argc - 5) && end > 0) {
    end = atom_getint(argv + argc - 1);
    argc--;
  }

  if(argc > 4 && !trim_check_number(argv + argc - 5)) {
    norm_mode = atom_getint(argv + argc - 1);
    argc--;
  }

  if(argc > 4 && !trim_check_number(argv + argc - 5)) {
    fade_type = atom_getint(argv + argc - 1);
    argc--;
  }
  
  if(argc > 4 && !trim_check_number(argv + argc - 5)) {
    normalize = atom_getfloat(argv + argc - 1);
    argc--;
  }
  
  if(argc > 4 && !trim_check_number(argv + argc - 5)) {
    pad_out_time = atom_getfloat(argv + argc - 1);
    argc--;
  }

  if(argc > 4 && !trim_check_number(argv + argc - 5)) {
    pad_in_time = atom_getfloat(argv + argc - 1);
    argc--;
  }

  if(sym == gensym("trim") && argc < 6) {
    pd_error(x, "not enough arguments to message %s", sym->s_name);
    return;
  }

  if(trim_check_db(argv + argc - 4, &in_db) || trim_check_db(argv + argc - 3, &out_db)) {
    pd_error(x, "expected number of 'off' to message %s", sym->s_name);
    return;
  }

  if(trim_check_number(argv + argc - 2) || trim_check_number(argv + argc - 1)) {
    pd_error(x, "expected number argument to message %s", sym->s_name);
    return;
  }

  fade_in_time = atom_getfloat(argv + argc - 2);
  fade_out_time = atom_getfloat(argv + argc - 1);
  argc -= 4;

  argc_halved = argc >> 1;

  num_buffers = buffer_multiple_names(in_buffer_names, out_buffer_names, lengths, argc, argv, in_place, &overall_length, &max_length);

  if(!num_buffers)
    return;

  fade_in = (t_int)(fade_in_time * sample_rate / 1000.);
  fade_out = (t_int)(fade_out_time * sample_rate / 1000.);
  pad_in = (t_int)(pad_in_time * sample_rate / 1000.);
  pad_out = (t_int)(pad_out_time * sample_rate / 1000.);

  fade_in = fade_in < 0 ? 0 : fade_in;
  fade_out = fade_out < 0 ? 0 : fade_out;
  pad_in = pad_in < 0 ? 0 : pad_in;
  pad_out = pad_out < 0 ? 0 : pad_out;

  samples[0] = (t_sample *)aligned_getbytes(sizeof(t_sample) * overall_length);
  temp_buf = (t_sample *)aligned_getbytes(sizeof(t_sample) * (max_length + pad_in + pad_out));

  if(!samples[0] || !temp_buf) {
    pd_error(x, "could not allocate temporary memory");
    aligned_freebytes(samples[0], sizeof(t_sample) * overall_length);
    aligned_freebytes(temp_buf, sizeof(t_sample) * (max_length + pad_in + pad_out));
  }

  for(i=0; i<num_buffers; i++) {
    samples[i] = samples[0] + offset;
    buffer_read(in_buffer_names[i], temp_buf, lengths[i]);

    for(j=0; j<lengths[i]; j++)
      samples[i][j] = temp_buf[j];
    
    offset += lengths[i];
  }

  norm_factor = calculate_trim(samples, lengths, max_length, overall_length,
			       num_buffers, in_db, out_db, sample_rate,
			       &trim_offset, &trim_length, normalize,
			       norm_mode);
  
  /* num_buffers = buffer_multiple_names(out_buffer_names, lengths_out, */
  /* 				      argc_halved, argv+argc_halved, */
  /* 				      &overall_length, &max_length); */

  /* if(!num_buffers) */
  /*   return; */
  
  if(trim_length==0)
    pd_error(x, "requested start trim level / end trim level never reached");
  else {
    for(i=0; i<num_buffers; i++)
      if(trim_write_buffer(out_buffer_names[i], samples[i], temp_buf,
			   trim_offset, trim_length, fade_in, fade_out, pad_in,
			   pad_out, norm_factor, lengths[i], sample_rate,
			   (t_fade_type)fade_type, end)) {}
	/* overall_error = true; */
  }

  aligned_freebytes(samples[0], sizeof(t_sample) * overall_length);
  aligned_freebytes(temp_buf, sizeof(t_sample) * (max_length + pad_in + pad_out));
  
  outlet_float(x->floatout, IRMANIP_TRIM);
}

static void irmanip_phase(t_irmanip *x, t_symbol *sym,
				  long argc, t_atom *argv)
{
  FFT_Setup *fft_setup;

  FFT_Split spectrum_1;
  FFT_Split spectrum_2;
  FFT_Split spectrum_3;

  t_float *in;
  /* t_float *filter_in; */
  /* non serve attualmente */
  t_sample *out_buf;

  t_symbol *filter = filter_retriever(x->deconvolve_filter_specifier);
  t_symbol *target = atom_getsymbol(argv++);
  t_symbol *source = atom_getsymbol(argv++);

  t_float filter_specifier[PIRO_MAX_SPECIFIER_ITEMS];
  t_float range_specifier[PIRO_MAX_SPECIFIER_ITEMS];

  t_float phase = 0.;
  t_float time_mul = 1.;
  t_float sample_rate = x->sample_rate;
  /* t_float deconvolve_phase; */

  /* t_float deconvdelay; */

  t_int mode = 0; /* */

  t_uint fft_size;
  t_uint fft_size_log2;
  t_uint i;

  /* long deconvolve_mode; */

  t_int source_length = buffer_length(source);
  t_int filter_length = buffer_length(filter);
  t_int max_length = source_length;

  fft_size = calculate_fft_size((long)(max_length * time_mul), &fft_size_log2);

  /* deconvolve_mode = x->deconvolve_mode; */
  /* deconvolve_phase = phase_retriever(x->deconvolve_phase); */
  /* deconvdelay = delay_retriever(&x->deconvolve_delay, fft_size, sample_rate); */

  fft_setup = create_setup(fft_size_log2);
  spectrum_1.realp = (t_sample *)aligned_getbytes(fft_size * 3 * sizeof(t_sample));
  spectrum_1.imagp = spectrum_1.realp + fft_size;
  spectrum_2.realp = spectrum_1.imagp + fft_size;

  /* filter_in = filter_length ? (t_float *)aligned_getbytes(filter_length * */
  /* 						  sizeof(t_float)) : NULL; */

  out_buf = spectrum_2.realp;
  in = (t_sample *)out_buf;

  /* if(!fft_setup || !spectrum_1.realp || (filter_length && !filter_in)) { */
  if(!fft_setup || !spectrum_1.realp) {
    pd_error(x, "could not allocate temporary memory");
    destroy_setup(fft_setup);
    aligned_freebytes(spectrum_1.realp, fft_size*3*sizeof(t_sample));
    /* aligned_freebytes(filter_in, filter_length ? filter_length*sizeof(t_float) : 0); */

    outlet_float(x->floatout, IRMANIP_ERROR);

    return;
  }

  buffer_read(source, in, fft_size);
  time_to_spectrum(fft_setup, in, source_length, spectrum_1, fft_size);
  power_spectrum(spectrum_1, fft_size, SPECTRUM_FULL);
  variable_phase_from_power_spectrum(fft_setup, spectrum_1, fft_size, phase, false);

  spectrum_to_time(fft_setup, out_buf, spectrum_1, fft_size, SPECTRUM_FULL);
  buffer_write(target, out_buf, fft_size, 1, 1.);

  destroy_setup(fft_setup);
  aligned_freebytes(spectrum_1.realp, fft_size*3*sizeof(t_sample));
  /* aligned_freebytes(filter_in, filter_length ? filter_length*sizeof(t_float) : 0); */
  
    outlet_float(x->floatout, IRMANIP_PHASE);
}

// MIMO (out-of-place only)

long irmanip_matrix_mimo(t_irmanip *x, t_matrix_complex *out, t_matrix_complex *in, t_sample regularization)
{
  t_uint m_dim = in->m_dim;
  t_uint n_dim = in->n_dim;
  t_uint i;
	
  MATRIX_REF_COMPLEX(out);
  MATRIX_REF_COMPLEX(temp_matrix_2);
	
  t_matrix_complex *temp_matrix_1;
  t_matrix_complex *temp_matrix_2;
  
  /* Check Dimensions / Size Output / Allocate Temporary Matrices */
  
  temp_matrix_1 = matrix_alloc_complex(n_dim, m_dim);
  temp_matrix_2 = matrix_alloc_complex(n_dim, n_dim);
	
  if (!temp_matrix_1 || !temp_matrix_2) {
    pd_error(x, "could not allocate matrix storage");
    return 1;
  }
	
  /* This call should not ever fail */
		
  if (matrix_new_size_complex(out, n_dim, n_dim)) {
    pd_error(x, "could not allocate matrix storage");
    return 1;
  }
	
  /* Dereference */
	
  MATRIX_DEREF(out);
  MATRIX_DEREF(temp_matrix_2);
	
  /* Calculate */
	
  matrix_conjugate_transpose_complex(temp_matrix_1, in);
  matrix_multiply_complex(out, temp_matrix_1, in);
	
  for (i = 0; i < n_dim; i++)
    MATRIX_ELEMENT(out, i, i) = CADD(MATRIX_ELEMENT(out, i, i), CSET(regularization, 0));
	
  matrix_choelsky_decompose_complex(temp_matrix_2, out);
  matrix_choelsky_solve_complex(out, temp_matrix_2, temp_matrix_1);
	
  /* Free */
	
  matrix_destroy_complex(temp_matrix_1);
  matrix_destroy_complex(temp_matrix_2);
	
  return 0;
}

long irmanip_mimo_deconvolution(t_irmanip *x, FFT_Split *impulses, t_uint fft_size, long sources, long receivers, t_sample *regularization)
{
  t_matrix_complex *in;
  t_matrix_complex *out;
	
  t_uint fft_size_halved = fft_size >> 1;
  t_uint i;
	
  long j, k;
	
  t_complex_sample out_val;
	
  MATRIX_REF_COMPLEX(in);
  MATRIX_REF_COMPLEX(out);
	
  in = matrix_alloc_complex(receivers, sources);
  out = matrix_alloc_complex(sources, receivers);
	
  MATRIX_DEREF(in);
  MATRIX_DEREF(out);
	
  /* Do DC */
	
  for(j = 0; j < receivers; j++)
    for (k = 0; k < sources; k++)
      MATRIX_ELEMENT(in, j, k) = CSET(impulses[j * sources + k].realp[0], 0);
	
  if(irmanip_matrix_mimo(x, out, in, regularization[0]))
    goto mimo_bail;
	
  for (j = 0; j < receivers; j++)
    for (k = 0; k < sources; k++)
      impulses[j * sources + k].realp[0] = CREAL(MATRIX_ELEMENT(out, j, k));
	
  /* Do Nyquist */
	
  for (j = 0; j < receivers; j++)
    for (k = 0; k < sources; k++)
      MATRIX_ELEMENT(in, j, k) = CSET(impulses[j * sources + k].imagp[0], 0);
	
  if(irmanip_matrix_mimo(x, out, in, regularization[fft_size_halved]))
    goto mimo_bail;
	
  for (j = 0; j < receivers; j++)
    for (k = 0; k < sources; k++)
      impulses[j * sources + k].imagp[0] = CREAL(MATRIX_ELEMENT(out, j, k));
	
  /* Do Other Bins */
	
  for (i = 1; i < fft_size_halved; i++) {		
    for (j = 0; j < receivers; j++) {
      for (k = 0; k < sources; k++)
	MATRIX_ELEMENT(in, j, k) = CSET(impulses[j * sources + k].realp[i], impulses[j * sources + k].imagp[i]);
    }
		
    if(irmanip_matrix_mimo(x, out, in, regularization[i]))
      break;
		
    for (j = 0; j < receivers; j++) {
      for (k = 0; k < sources; k++) {
	out_val = MATRIX_ELEMENT(out, j, k);
	impulses[j * sources + k].realp[i] = CREAL(out_val);
	impulses[j * sources + k].imagp[i] = CIMAG(out_val);
      }
    }
  }
	
  if(i == fft_size_halved) {
    matrix_destroy_complex(in);
    matrix_destroy_complex(out);
	
    return 0;
  }
	
mimo_bail:
	
  matrix_destroy_complex(in);
  matrix_destroy_complex(out);
	
  return 1;
	
}

void irmanip_mimo(t_irmanip *x, t_symbol *sym, long argc,
		       t_atom *argv)
{
  FFT_Setup *fft_setup;

  FFT_Split impulses[128];

  t_symbol *in_buffer_names[128];
  t_symbol *out_buffer_names[128];

  t_float filter_specifier[PIRO_MAX_SPECIFIER_ITEMS];

  t_sample *temp_buffer;
  t_sample *regularization;
  /* float *temp_buffer; */

  t_int lengths[128];
  /* t_int lengths_out[128]; */
  /* optional */

  double sample_rate = sys_getsr();
  double time_mul = 1.;
  double deconvolve_delay;

  long receivers;
  t_int sources;

  t_uint fft_size;
  t_uint fft_size_log2;

  t_int overall_length = 0;
  t_int max_length = 0;
  t_int length;
  t_int num_buffers = 0;
  /* t_int num_buffers_out = 0; */
  t_int i, j;

  int error;
  t_int in_place = 1;
  /* long read_chan = x->read_chan - 1;  */
  long write_chan = x->write_chan - 1;

  t_bool overall_error = false;

  /* int argc_halved = argc; */

  if(sym == gensym("mimoto")) {
    in_place = 0;
    /* argc_halved = argc >> 1; */
  }

  /* Get number of sources / receivers */

  if(argc < 2) {
    pd_error(x, "not enough arguments to message %s", sym->s_name);
    return;
  }

  sources = atom_getint(argv++);
  argc--;

  if (sources <= 0) {
    pd_error(x, "not enough arguments to message %s", sym->s_name);
    return;
  }

  if(argc && argv->a_type == A_FLOAT) {
    time_mul = atom_getfloat(argv++);
    argc--;

    if(time_mul < 1.) {
      post("time multiplier cannot be less than 1 (using 1)");
      time_mul = 1.;
    }
  }

  receivers = in_place ? argc / sources : argc / (sources * 2);
  num_buffers = in_place ? argc : argc / 2;
  /* num_buffers_out = num_buffers; */

  if(sources * receivers != num_buffers) {
    pd_error(x, "number of specified buffers is not divisible by the number of sources - additional or missing IRs");
    return;
  }

  /* Check buffers, storing names and lengths +  calculate total / largest length */

  num_buffers = buffer_multiple_names(in_buffer_names, out_buffer_names, lengths, argc, argv, in_place, &overall_length, &max_length);

  /* post("num_buffers = %d argc = %d", num_buffers, argc); */

  if(!num_buffers)
    return;

  /* for(i=0; i<num_buffers; i++) */
  /*   post("in_buffer = %s out_buffer = %s", in_buffer_names[i]->s_name, */
  /* 	 out_buffer_names[i]->s_name); */

  /* Calculate fft size */

  fft_size = calculate_fft_size((t_uint) (max_length * time_mul),
  				&fft_size_log2);

  
  deconvolve_delay = delay_retriever(&x->deconvolve_delay, fft_size,
  				     sample_rate);

  /* Check length of buffers for writing */

  if(!x->resize) {
    for(i = 0; i<num_buffers; i++) {
      if(buffer_length(out_buffer_names[i]) < (t_int) fft_size) {
  	pd_error(x, "buffer %s is not long enough to complete write (no buffers altered)", out_buffer_names[i]->s_name);
  	return;
      }
    }
  }

  /* Allocate Resources */

  fft_setup = create_setup(fft_size_log2);
  temp_buffer = (t_sample *)aligned_getbytes(sizeof(t_sample) * fft_size * 2);
  regularization = temp_buffer + fft_size;
  impulses[0].realp = (t_sample *)aligned_getbytes(sizeof(t_sample) * fft_size * sources * receivers);

  /* Check Memory Allocations */

  if(!fft_setup || !temp_buffer || !impulses[0].realp) {
    pd_error(x, "could not allocate temporary memory for processing");

    destroy_setup(fft_setup);
    aligned_freebytes(impulses[0].realp, sizeof(t_sample) * fft_size * sources * receivers);
    aligned_freebytes(temp_buffer, sizeof(t_sample) * fft_size * 2);
    
    outlet_float(x->floatout, IRMANIP_ERROR);
    
    return;
  }

  /* Set pointers for impulses */

  for(i = 0; i < sources * receivers; i++) {
    impulses[i].realp = impulses[0].realp + (fft_size * i);
    impulses[i].imagp = impulses[i].realp + (fft_size >> 1);
  }

  /* Do transforms in */

  for(i = 0; i < receivers * sources; i++) {
    length = buffer_read(in_buffer_names[i], temp_buffer, fft_size);
    time_to_halfspectrum(fft_setup, temp_buffer, length, impulses[i], fft_size);
  }

  /* Prepare regularisation */

  fill_power_array_specifier(filter_specifier, x->deconvolve_filter_specifier,
			     x->deconvolve_num_filter_specifiers);
  make_freq_dependent_power_array(regularization, filter_specifier, fft_size,
				  sample_rate, 0);
  

  /* Deconvolve */

  if(!irmanip_mimo_deconvolution(x, impulses, fft_size, sources, receivers, regularization)) {
    
    /* Do transforms out */

    for(i = 0; i < receivers * sources; i++) {
      delay_spectrum(impulses[i], fft_size, SPECTRUM_REAL, deconvolve_delay);
      spectrum_to_time(fft_setup, temp_buffer, impulses[i], fft_size,
      		       SPECTRUM_REAL);
      buffer_write(out_buffer_names[i], temp_buffer, fft_size, 1, 1.);

      
      /* buffer_write_error((t_object *) x, out_buffer_names[i], error); */

      /* if(error) */
      /* 	overall_error = true; */
    }
  }

  destroy_setup(fft_setup);
  aligned_freebytes(impulses[0].realp, sizeof(t_sample) * fft_size * sources * receivers);
  aligned_freebytes(temp_buffer, sizeof(t_sample) * fft_size * 2);

  outlet_float(x->floatout, IRMANIP_INVERT);
}

static void irmanip_invert(t_irmanip *x, t_symbol *sym,
				long argc, t_atom *argv)
{
  FFT_Setup *fft_setup;

  FFT_Split spectrum_1;
  FFT_Split spectrum_2;
  FFT_Split spectrum_3;

  t_sample *out_buf;
  t_sample *in_temp;
  t_float *filter_in;

  t_symbol *target = atom_getsymbol(argv++);
  t_symbol *source = atom_getsymbol(argv++);
  
  t_symbol *filter = filter_retriever(x->deconvolve_filter_specifier);
  
  int deconvfilter_length = 0;
  t_garray *buf;
  t_word *buf_samples;

  t_float filter_specifier[PIRO_MAX_SPECIFIER_ITEMS];
  t_float range_specifier[PIRO_MAX_SPECIFIER_ITEMS];

  t_float time_mul = 1.;
  
  if(argc>2)
    time_mul = atom_getfloat(argv++);
  
  t_float sample_rate = x->sample_rate;
  t_float deconvolve_phase = phase_retriever(x->deconvolve_phase);
  
  t_float deconvdelay;

  t_int filter_length = buffer_length(filter);
  t_int source_length = buffer_length(source);

  t_uint fft_size;
  t_uint fft_size_log2;

  long deconvolve_mode = x->deconvolve_mode;

  fft_size = calculate_fft_size((t_int)(source_length * time_mul),
				&fft_size_log2);
  fft_setup = create_setup(fft_size_log2);

  deconvdelay = delay_retriever(&x->deconvolve_delay, fft_size, sample_rate);

  spectrum_1.realp = (t_sample *)aligned_getbytes(fft_size * 4 * sizeof(t_sample));
  spectrum_1.imagp = spectrum_1.realp + (fft_size >> 1);
  spectrum_2.realp = spectrum_1.imagp + (fft_size >> 1);
  spectrum_2.imagp = spectrum_2.realp + (fft_size >> 1);
  spectrum_3.realp = spectrum_2.imagp + (fft_size >> 1);
  spectrum_3.imagp = spectrum_3.realp + fft_size;

  attach_array(filter, &buf, &buf_samples, &deconvfilter_length);
  
  filter_in = deconvfilter_length > 0 ?
    (t_float *)aligned_getbytes(fft_size * sizeof(t_float)) : NULL;

  in_temp = (t_sample *)spectrum_3.realp;
  out_buf = spectrum_3.realp;

  if(!fft_setup || !spectrum_1.realp || (deconvfilter_length > 0 && !filter_in)) {
    pd_error(x, "could not allocate temporary memory");
    destroy_setup(fft_setup);
    aligned_freebytes(spectrum_1.realp, fft_size*4*sizeof(t_sample));
    aligned_freebytes(filter_in, deconvfilter_length > 0 ?
	      fft_size*sizeof(t_float) : 0);

    outlet_float(x->floatout, IRMANIP_ERROR);
    
    return;
  }

  spike_spectrum(spectrum_1, fft_size, SPECTRUM_REAL, deconvdelay);
  buffer_read(source, in_temp, fft_size);
  time_to_halfspectrum(fft_setup, in_temp, fft_size, spectrum_2, fft_size);
  
  fill_power_array_specifier(filter_specifier, x->deconvolve_filter_specifier,
  			     x->deconvolve_num_filter_specifiers);
  fill_power_array_specifier(range_specifier, x->deconvolve_range_specifier,
  			     x->deconvolve_num_range_specifiers);
  buffer_read(filter, (t_sample *)filter_in, fft_size);
  deconvolve(fft_setup, spectrum_1, spectrum_2, spectrum_3, filter_specifier, range_specifier, 0, filter_in, fft_size, fft_size, SPECTRUM_REAL, (t_filter_type)deconvolve_mode, deconvolve_phase, 0., sample_rate);

  spectrum_to_time(fft_setup, out_buf, spectrum_1, fft_size, SPECTRUM_REAL);
  buffer_write(target, out_buf, fft_size, 1, 1.);

  destroy_setup(fft_setup);
  aligned_freebytes(spectrum_1.realp, fft_size*4*sizeof(t_sample));
  aligned_freebytes(filter_in, deconvfilter_length ? fft_size*sizeof(t_float) : 0);
  
  outlet_float(x->floatout, IRMANIP_INVERT);
}

static void irmanip_average(t_irmanip *x, t_symbol *sym,
				 long argc, t_atom *argv)
{
  FFT_Setup *fft_setup;

  FFT_Split spectrum_1;
  FFT_Split spectrum_2;

  t_sample *temp;

  t_symbol *target;
  t_symbol *buffer_names[128];

  t_float sample_rate = 0.;

  /* smooth attributes */
  long smooth_mode = 0;
  long num_smooth = 0;
  double smooth[2] = {0, 0.9999};

  t_int lengths[128];

  t_uint fft_size;
  t_uint fft_size_log2;

  t_int num_buffers = 0;
  t_int max_length;
  t_int read_length;
  t_int overall_length;
  t_int i, j;

  t_int err = 0;

  if(!argc) {
    pd_error(x, "not enough arguments to message %s", sym->s_name);
  }

  target = atom_getsymbol(argv++);
  argc--;

  if(argc && argv->a_type == A_FLOAT) {
    smooth_mode = (t_int)atom_getint(argv++);
    argc--;

    if(smooth_mode < 0) {
      pd_error(x, "smooth set to zero");
      smooth_mode = 0;
    }

    if(smooth_mode > 3) {
      pd_error(x, "smooth set to zero");
      smooth_mode = 0;
    }
  }

  if(argc && argv->a_type == A_FLOAT) {
    smooth[0] = atom_getfloat(argv++);
    argc--;
    num_smooth++;
  }

  if(argc && argv->a_type == A_FLOAT) {
    smooth[1] = atom_getfloat(argv++);
    argc--;
    num_smooth++;
  }
  
  num_buffers = buffer_multiple_names(buffer_names, NULL, lengths, argc, argv, 0, &overall_length, &max_length);
  
  if(!num_buffers)
    return;

  fft_size = calculate_fft_size((t_uint)max_length,&fft_size_log2);

  fft_setup = create_setup(fft_size_log2);

  temp = (t_sample *)aligned_getbytes(fft_size * 5 * sizeof(t_sample));
  spectrum_1.realp = temp + fft_size;
  spectrum_1.imagp = spectrum_1.realp + fft_size;
  spectrum_2.realp = spectrum_1.imagp + fft_size;
  spectrum_2.imagp = spectrum_2.realp + fft_size;

  if(!fft_setup || !temp) {
    pd_error(x, "could not allocate temporary memory");
    destroy_setup(fft_setup);
    aligned_freebytes(temp, fft_size*5*sizeof(t_sample));
    
    outlet_float(x->floatout, IRMANIP_ERROR);
  }

  for(j=0; j<(t_int)fft_size; j++) {
    spectrum_1.realp[j] = 0.;
    spectrum_1.imagp[j] = 0.;
  }

  for(i=0; i<num_buffers; i++) {
    read_length = buffer_read(buffer_names[i], temp, fft_size); /* */
    time_to_spectrum(fft_setup, temp, read_length, spectrum_2, fft_size);
    power_spectrum(spectrum_2, fft_size, SPECTRUM_FULL);

    for(j=0; j<(t_int)fft_size; j++)
      spectrum_1.realp[j] += spectrum_2.realp[j] / (t_sample)num_buffers;
  }
  
  if(num_smooth) {
    smooth_power_spectrum(spectrum_1, smooth_mode, fft_size, num_smooth > 1 ? smooth[0] : 0., num_smooth > 1 ? smooth[1] : smooth[0]);
  }

  variable_phase_from_power_spectrum(fft_setup, spectrum_1, fft_size, phase_retriever(x->deconvolve_phase), false); /* fase? */
  spectrum_to_time(fft_setup, temp, spectrum_1, fft_size, SPECTRUM_FULL);
  err = buffer_write(target, temp, fft_size, 1, 1.);
  if(!err)
    pd_error(x, "%s array error", target->s_name);

  
  destroy_setup(fft_setup);
  aligned_freebytes(temp, fft_size*5*sizeof(t_sample));
  
  outlet_float(x->floatout, IRMANIP_AVERAGE);
}

static void irmanip_free(t_irmanip *x)
{
  aligned_freebytes(x->deconvolve_filter_specifier,
	    PIRO_MAX_SPECIFIER_ITEMS*sizeof(t_atom));
  aligned_freebytes(x->deconvolve_range_specifier,
	    PIRO_MAX_SPECIFIER_ITEMS*sizeof(t_atom));
}

static void *irmanip_new(t_symbol *s, short argc, t_atom *argv)
{
  t_irmanip *x = (t_irmanip *)pd_new(irmanip_class);

  /* deconvolution attributes & methods */
  x->deconvolve_mode = 0; /* FILTER_REGULARISATION */
  x->deconvolve_filter_specifier = (t_atom *)aligned_getbytes(PIRO_MAX_SPECIFIER_ITEMS*sizeof(t_atom));
  x->deconvolve_range_specifier = (t_atom *)aligned_getbytes(PIRO_MAX_SPECIFIER_ITEMS*sizeof(t_atom)); // to check
  if(!x->deconvolve_filter_specifier || !x->deconvolve_range_specifier)
    pd_error(x,"could not allocate space for attribute storage");
  irmanip_deconvrange(x, 0, 0, 0);
  irmanip_deconvfilter(x, 0, 0, 0);
  irmanip_deconvphase(x, 0, 0, 0);
  irmanip_deconvdelay(x, 0, 0, 0);

  x->write_chan = 1;
  x->resize = 1;
  
  x->sample_rate = sys_getsr();
  if(!x->sample_rate)
    x->sample_rate = 44100.0;

  x->floatout = outlet_new(&x->x_obj, &s_float);

  return (void *)x;
}

void irmanip_setup(void)
{
#if PD_FLOATSIZE == 32
  irmanip_class = class_new(gensym("irmanip"),
				 (t_newmethod)irmanip_new,
				 (t_method)irmanip_free,
				 sizeof(t_irmanip),
				 CLASS_DEFAULT, 
				 A_GIMME, 0);
#elif PD_FLOATSIZE == 64
  irmanip_class = class_new64(gensym("irmanip"),
				   (t_newmethod)irmanip_new,
				   (t_method)irmanip_free,
				   sizeof(t_irmanip),
				   CLASS_DEFAULT, 
				   A_GIMME, 0);
#else
#error [irmanip]: invalid FLOATSIZE: must be 32 or 64
#endif
  
  class_addmethod(irmanip_class,
		  (t_method)irmanip_average, gensym("average"),
		  A_GIMME, 0);
  class_addmethod(irmanip_class,
		  (t_method)irmanip_phase, gensym("phase"),
		  A_GIMME, 0);
  
  
  class_addmethod(irmanip_class,
		  (t_method)irmanip_trim, gensym("trim"),
		  A_GIMME, 0);
  class_addmethod(irmanip_class,
		  (t_method)irmanip_trim, gensym("trimto"),
		  A_GIMME, 0);


  class_addmethod(irmanip_class,
		  (t_method)irmanip_invert, gensym("invert"),
		  A_GIMME, 0);

  class_addmethod(irmanip_class,
		  (t_method)irmanip_mimo, gensym("mimo"),
		  A_GIMME, 0);
  class_addmethod(irmanip_class,
		  (t_method)irmanip_mimo, gensym("mimoto"),
		  A_GIMME, 0);

  class_addmethod(irmanip_class,
		  (t_method)irmanip_deconvmode, gensym("deconvmode"),
		  A_FLOAT, 0);
  class_addmethod(irmanip_class,
		  (t_method)irmanip_deconvrange, gensym("deconvrange"),
		  A_GIMME, 0);

  class_addmethod(irmanip_class,
		  (t_method)irmanip_deconvfilter,
		  gensym("deconvfilter"),A_GIMME, 0);
           
  class_addmethod(irmanip_class,
		  (t_method)irmanip_deconvphase, gensym("deconvphase"),
		  A_GIMME, 0);

  class_addmethod(irmanip_class,
		  (t_method)irmanip_deconvdelay, gensym("deconvdelay"),
		  A_GIMME, 0);
}
