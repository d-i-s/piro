#include "z_core.h"

double db_to_pow(double db)
{
  return pow(10., db / 10.);
}

void db_to_pow_array(t_sample *in, t_uint length)
{
  t_uint i;

  for(i=0;i<length;i++)
    in[i] = db_to_pow(in[i]);
}

double pow_to_db(double pow)
{
  double db;

  if(!pow)
    return PIRO_DB_MIN;

  db = 10. * log10(pow);

  if(db < PIRO_DB_MIN)
    db = PIRO_DB_MIN;

  return db;
}

double db_to_a(double db)
{
  return pow(10., db / 20.);
}




void spike_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format, double spike)
{
  long double spike_const = (long double) (2. * M_PI) * (double) (fft_size - spike) / ((double) fft_size);
  long double phase;
  t_uint i;

  spectrum.realp[0] = 1;

  if(format == SPECTRUM_FULL)
    spectrum.imagp[0] = 0.;
  else
    spectrum.imagp[0] = 1.;

  for(i=1; i<fft_size>>1; i++) {
    phase = spike_const * i;
    spectrum.realp[i] = cosl(phase);
    spectrum.imagp[i] = sinl(phase);
  }

  if(format == SPECTRUM_FULL) {
    spectrum.realp[i] = 1.;
    spectrum.imagp[i] = 0.;

    for(i++; i<fft_size; i++) {
      spectrum.realp[i] = spectrum.realp[fft_size - i];
      spectrum.imagp[i] = -spectrum.imagp[fft_size - i];
    }
  }
}

void delay_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format, double delay)
{
  long double delay_const = (long double) (2.0 * M_PI) * (double) -delay /
    ((double) fft_size);
  long double phase;
  double a, b, c, d;
  t_uint i;

  if(!delay)
    return;
  
  for(i=1;i<fft_size>>1;i++) {
    phase = delay_const * i;

    a = spectrum.realp[i];
    b = spectrum.imagp[i];
    c = cosl(phase);
    d = sinl(phase);

    spectrum.realp[i] = (a * c) - (b * d);
    spectrum.imagp[i] = (b * c) - (a * d);
  }

  if(format==SPECTRUM_FULL) {
    for(i++; i<fft_size; i++) {
      spectrum.realp[i] = spectrum.realp[fft_size - i];
      spectrum.imagp[i] = -spectrum.imagp[fft_size - i];
    }
  }
}

void make_freq_dependent_power_array(t_sample *power_array, t_sample *specifier_array, t_uint fft_size, t_sample sample_rate, t_sample db_offset)
{
  t_uint fft_size_halved = fft_size >> 1;
  t_uint num_items = 0;
  t_uint list_pos;
  t_uint i;

  t_sample freq_mul = sample_rate / ((t_sample) fft_size);
  double prev_log_freq = -HUGE_VAL;
  double next_log_freq;
  t_sample bin_log_freq = 0.;
  double gradient;
  double offset;

  for(i=0; i<PIRO_MAX_SPECIFIER_ITEMS; i++) {
    if(isinf(specifier_array[i])) 
      break;
  }

  num_items = i;

  if(num_items<=2) {
    t_sample pow_val = num_items == 1 ?
      db_to_pow(specifier_array[0] + db_offset)
      : db_to_pow(specifier_array[1] + db_offset);
    
    for(i=0; i<fft_size; i++)
      power_array[i] = pow_val;

    return;
  }

  num_items >>= 1;
  power_array[0] = specifier_array[1] + db_offset;

  for(i=1, list_pos = 0, gradient = 0., offset = specifier_array[1],
  	next_log_freq = log(specifier_array[0]); i<fft_size_halved+1; i++) {
    
    bin_log_freq = log(i * freq_mul);

    if(bin_log_freq > next_log_freq) {
      for(; bin_log_freq>next_log_freq && list_pos < num_items; list_pos++,
    	    prev_log_freq = next_log_freq,
    	    next_log_freq = log(specifier_array[list_pos << 1]));

      if(list_pos == num_items) {
    	gradient = 0.;
    	offset = specifier_array[(list_pos << 1) - 1];
    	next_log_freq = HUGE_VAL;
      }
      else {
    	gradient = (specifier_array[(list_pos << 1) + 1] -
    		    specifier_array[(list_pos << 1) - 1]) /
    	  (next_log_freq - prev_log_freq);
    	offset = specifier_array[(list_pos << 1) - 1] -
    	  (prev_log_freq * gradient);
      }
    }
    
    power_array[i] = bin_log_freq * gradient + offset + db_offset;
  }
  
  db_to_pow_array(power_array, fft_size_halved+1);

  for(i=fft_size_halved+1; i<fft_size; i++)
    power_array[i] = power_array[fft_size - i];
}

void deconvolve_regularised_zero_phase(FFT_Split spectrum_1, FFT_Split spectrum_2, t_float *beta_in, t_uint fft_size, t_spectrum_format format)
{
  t_sample *real1 = spectrum_1.realp;
  t_sample *imag1 = spectrum_1.imagp;
  t_sample *real2 = spectrum_2.realp;
  t_sample *imag2 = spectrum_2.imagp;

  t_float a, b, c, d, e;

  t_uint from = format == SPECTRUM_REAL ? 1 : 0;
  t_uint to = format == SPECTRUM_REAL ? fft_size >> 1 : fft_size;
  t_uint i;

  if(format == SPECTRUM_REAL) {
    a = real1[0];
    b = imag1[0];
    c = real2[0];
    d = imag2[0];

    real1[0] = a*c / (c*c + beta_in[0]);
    imag1[0] = b*d / (d*d + beta_in[fft_size>>1]);
  }

  for(i=from; i<to; i++) {
    a = real1[i];
    b = imag1[i];
    c = real2[i];
    d = imag2[i];

    e = 1. / (((c*c) + (d*d)) + beta_in[i]);

    real1[i] = ((a*c) + (b*d)) * e;
    imag1[i] = ((b*c) - (a*d)) * e;
  }
}

void deconvolve_with_filter(FFT_Split spectrum_1, FFT_Split spectrum_2, FFT_Split filter_spectrum, t_uint fft_size, t_spectrum_format format)
{
  t_sample *real1 = spectrum_1.realp;
  t_sample *imag1 = spectrum_1.imagp;
  t_sample *real2 = spectrum_2.realp;
  t_sample *imag2 = spectrum_2.imagp;
  t_sample *real3 = filter_spectrum.realp;
  t_sample *imag3 = filter_spectrum.imagp;

  t_complex_sample c1, c2, c3;

  t_float mul;

  t_uint from = format == SPECTRUM_REAL ? 1 : 0;
  t_uint to = format == SPECTRUM_REAL ? fft_size >> 1 : fft_size;
  t_uint i;

  if(format == SPECTRUM_REAL) {
    real1[0] = (real1[0] / real2[0]) * real3[0];
    imag1[0] = (imag1[0] / imag2[0]) * imag3[0];

    real1[0] = isinf(real1[0]) ? 0. : real1[0];
    imag1[0] = isinf(imag1[0]) ? 0. : imag1[0];
  }

  for(i=from; i<to; i++) {
    c1 = CSET(real1[i], imag1[i]);
    c2 = CSET(real2[i], imag2[i]);
    c3 = CSET(real3[i], imag3[i]);

    c1 = CMUL(CMUL(c1, CONJ(c2)), c3);

    mul = 1. / CABS_SQ(c2);
    mul = isinf(mul) ? 0. : mul;

    real1[i] = CREAL(c1) * mul;
    imag1[i] = CIMAG(c1) * mul;
  }
}

void deconvolve_with_amp_filter(FFT_Split spectrum_1, FFT_Split spectrum_2, t_float *filter_amps, t_uint fft_size, t_spectrum_format format)
{
  t_sample *real1 = spectrum_1.realp;
  t_sample *imag1 = spectrum_1.imagp;
  t_sample *real2 = spectrum_2.realp;
  t_sample *imag2 = spectrum_2.imagp;

  t_complex_sample c1, c2;

  t_sample mul;

  t_uint from = format == SPECTRUM_REAL ? 1 : 0;
  t_uint to = format == SPECTRUM_REAL ? fft_size >> 1 : fft_size;
  t_uint i;

  if(format == SPECTRUM_REAL) {
    real1[0] = (real1[0] / real2[0]) * filter_amps[0];
    imag1[0] = (imag1[0] / imag2[0]) * filter_amps[fft_size >> 1];

    real1[0] = isinf(real1[0]) ? 0. : real1[0];
    imag1[0] = isinf(imag1[0]) ? 0. : imag1[0];
  }

  for(i=from; i<to; i++) {
    c1 = CSET(real1[i], imag1[i]);
    c2 = CSET(real2[i], imag2[i]);

    c1 = CMUL(c1, CONJ(c2));

    mul = filter_amps[i] / CABS_SQ(c2);
    mul = isinf(mul) ? 0. : mul;

    real1[i] = CREAL(c1) * mul;
    imag1[i] = CIMAG(c1) * mul;
  }
}

void deconvolve_clip_zero_phase(FFT_Split spectrum_1, FFT_Split spectrum_2, t_float *clip_min, t_float *clip_max, t_uint fft_size, t_spectrum_format format)
{
  t_sample *real1 = spectrum_1.realp;
  t_sample *real2 = spectrum_2.realp;
  t_sample *imag1 = spectrum_1.imagp;
  t_sample *imag2 = spectrum_2.imagp;

  t_sample a, b, c, d, e, f = 0., min_sq, max_sq;

  t_uint from = format == SPECTRUM_REAL ? 1 : 0;
  t_uint to = format == SPECTRUM_REAL ? fft_size >> 1 : fft_size;
  t_uint i;

  if(format == SPECTRUM_REAL) {
    a = real1[0];
    b = imag1[0];
    c = real2[0];
    d = imag2[0];

    min_sq = clip_min[0];
    max_sq = clip_max[0];
    c = (c * c < min_sq) ? sqrt(min_sq) : c;
    c = (c * c > max_sq) ? sqrt(max_sq) : c;

    min_sq = clip_min[fft_size >> 1];
    max_sq = clip_max[fft_size >> 1];
    d = (d * d < min_sq) ? sqrt(min_sq) : d;
    d = (d * d > max_sq) ? sqrt(max_sq) : d;

    real1[0] = a / c;
    imag1[0] = b / d;
  }

  for(i=from; i<to; i++) {
    min_sq = clip_min[i];
    max_sq = clip_max[i];

    a = real1[i];
    b = imag1[i];
    c = real2[i];
    d = imag2[i];

    e = ((c * c) + (d * d));

    f = e < min_sq ? min_sq : e;
    f = f > max_sq ? max_sq : f;

    if(f != e) {
      e = sqrt(f / e);
      c *= e;
      d *= e;
      e = ((c * c) + (d * d));
    }

    e = 1. / e;

    real1[i] = ((a * c) + (b * d)) * e;
    imag1[i] = ((b * c) - (a * d)) * e;
    
  }
}

void deconvolve_zero_phase(FFT_Split spectrum_1, FFT_Split spectrum_2, FFT_Split filter_spectrum, t_float *filter_specifier, t_float *range_specifier, t_float filter_db_offset, t_uint fft_size, t_spectrum_format format, t_filter_type mode, t_float sample_rate)
{
  t_uint i;

  make_freq_dependent_power_array(filter_spectrum.realp, filter_specifier, fft_size, sample_rate, filter_db_offset);

  switch(mode)
    {
      
    case FILTER_REGULARISATION:
      for(i=0; i<fft_size; i++)
	filter_spectrum.imagp[i] = 0.;
      
      deconvolve_regularised_zero_phase(spectrum_1, spectrum_2, filter_spectrum.realp, fft_size, format);
      break;
      
    case FILTER_CLIP:
      make_freq_dependent_power_array(filter_spectrum.imagp, range_specifier, fft_size, sample_rate, 0.);

      for(i=0; i<(fft_size>>1); i++)
	filter_spectrum.imagp[i] *= filter_spectrum.realp[i];

      deconvolve_clip_zero_phase(spectrum_1, spectrum_2, filter_spectrum.realp, filter_spectrum.imagp, fft_size, format);
      break;

    case FILTER_FILTER:
      for(i=0; i<fft_size; i++)
	filter_spectrum.realp[i] = sqrt(filter_spectrum.realp[i]);

      deconvolve_with_amp_filter(spectrum_1, spectrum_2, filter_spectrum.realp, fft_size, format);
      break;
    }
}

void zero_phase_from_power_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format)
{
  t_uint to = format == SPECTRUM_REAL ? (fft_size >> 1) : (fft_size >> 1) +
    1;
  t_uint i;

  for(i=0;i<to;i++)
    spectrum.realp[i] = sqrt(spectrum.realp[i]);

  if(format==SPECTRUM_FULL) {
    for(i=to;i<fft_size;i++)
      spectrum.realp[i] = spectrum.realp[fft_size - i];
  }
}

void linear_phase_from_power_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format)
{
  zero_phase_from_power_spectrum(spectrum,fft_size,format);
  delay_spectrum(spectrum,fft_size,format,(double)(fft_size>>1));
}

void minimum_phase_components_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size)
{
  double scale = 1. / fft_size;
  double min_power = db_to_pow(-1000.);

  t_uint inexact;
  t_uint fft_size_halved = fft_size >> 1;
  t_uint fft_size_log2 = int_log2(fft_size, &inexact);
  t_uint i;

  for(i=0;i<fft_size_halved+1;i++) {
    spectrum.realp[i] = spectrum.realp[i] < min_power ? min_power :
      spectrum.realp[i];
    spectrum.realp[i] = 0.5 * log(spectrum.realp[i]);
  }

  for(i=fft_size_halved+1;i<fft_size;i++)
    spectrum.realp[i] = spectrum.realp[fft_size-i];

  do_ifft(&spectrum, fft_setup, fft_size_log2);

  for(i=1;i<fft_size_halved;i++) {
    spectrum.realp[i] += spectrum.realp[i];
    spectrum.imagp[i] += spectrum.imagp[i];
  }

  for(i=fft_size_halved+1;i<fft_size;i++) {
    spectrum.realp[i] = 0.;
    spectrum.imagp[i] = 0.;
  }

  for(i=0;i<fft_size_halved+1;i++) {
    spectrum.realp[i] *= scale;
    spectrum.imagp[i] *= scale;
  }

  do_fft(&spectrum, fft_setup, fft_size_log2);
}

void minimum_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size)
{
  t_uint fft_size_halved = fft_size >> 1;
  t_uint i;

  minimum_phase_components_from_power_spectrum(fft_setup, spectrum, fft_size);

  for(i=0;i<fft_size_halved+1;i++) {
    t_complex_sample c1 = CSET(spectrum.realp[i], spectrum.imagp[i]);
    t_complex_sample c2 = CEXP(c1);

    spectrum.realp[i] = CREAL(c2);
    spectrum.imagp[i] = CIMAG(c2);
  }

  for(i=fft_size_halved+1;i<fft_size;i++) {
    spectrum.realp[i] = spectrum.realp[fft_size - i];
    spectrum.imagp[i] = -spectrum.imagp[fft_size - i];
  }
}

void noncausal_maximum_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size)
{
  t_uint i;

  minimum_phase_from_power_spectrum(fft_setup, spectrum, fft_size);

  for(i=0;i<fft_size;i++)
    spectrum.imagp[i] = -spectrum.imagp[i];
}

void maximum_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size)
{
  noncausal_maximum_phase_from_power_spectrum(fft_setup, spectrum, fft_size);
  delay_spectrum(spectrum, fft_size, SPECTRUM_FULL, -1.);
}

void mixed_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size, double phase, t_bool zero_center)
{
  t_uint fft_size_halved = fft_size >> 1;
  t_uint i;

  double linphase_mul;
  double minphase_mul;
  double interp_phase;
  double amp;

  phase = phase < 0 ? 0. : phase;
  phase = phase > 1 ? 1. : phase;

  minphase_mul = 1 - (2 * phase);
  linphase_mul = zero_center ? 0. : (phase <= 0.5 ? -(2 * M_PI * phase) :
				     (-2 * M_PI * (phase - 1. / fft_size)));

  minimum_phase_components_from_power_spectrum(fft_setup, spectrum, fft_size);

  for(i=0; i<fft_size_halved;i++) {
    amp = exp(spectrum.realp[i]);
    interp_phase = linphase_mul * i + minphase_mul * spectrum.imagp[i];
    spectrum.realp[i] = amp * cos(interp_phase);
    spectrum.imagp[i] = amp * sin(interp_phase);
  }

  spectrum.realp[fft_size_halved] = exp(spectrum.realp[fft_size_halved]);
  spectrum.imagp[fft_size_halved] = 0.;

  for(i=fft_size_halved+1;i<fft_size;i++) {
    spectrum.realp[i] = spectrum.realp[fft_size - i];
    spectrum.imagp[i] = -spectrum.imagp[fft_size - i];
  }
}

void variable_phase_from_power_spectrum(FFT_Setup *fft_setup, FFT_Split spectrum, t_uint fft_size, double phase, t_bool zero_center)
{
  if(phase==0.) {
    minimum_phase_from_power_spectrum(fft_setup, spectrum, fft_size);
    return;
  }

  if(phase==0.5) {
    if(zero_center)
      zero_phase_from_power_spectrum(spectrum,fft_size,SPECTRUM_FULL);
    else
      linear_phase_from_power_spectrum(spectrum,fft_size,SPECTRUM_FULL);
    return;
  }

  if(phase==1.) {
    if(zero_center)
      noncausal_maximum_phase_from_power_spectrum(fft_setup, spectrum, fft_size);
    else
      maximum_phase_from_power_spectrum(fft_setup, spectrum,fft_size);
    return;
  }

  mixed_phase_from_power_spectrum(fft_setup, spectrum, fft_size, phase,
				  zero_center);
}

void power_spectrum(FFT_Split spectrum, t_uint fft_size, t_spectrum_format format)
{
  t_sample *real = spectrum.realp;
  t_sample *imag = spectrum.imagp;
  t_uint to = format == SPECTRUM_REAL ? fft_size >> 1 : fft_size;
  t_uint i;

  if(format == SPECTRUM_REAL)
    spectrum.imagp[0] = 0.;

  for(i=0; i<to; i++) {
    real[i] = (real[i] * real[i]) + (imag[i] * imag[i]);
    imag[i] = 0.;
  }
}

void make_regularisation_filter(FFT_Setup *fft_setup, FFT_Split denominator_spectrum, FFT_Split filter_spectrum, t_sample *beta_in, t_uint fft_size, t_spectrum_format format, double phase)
{
  t_sample *real1 = denominator_spectrum.realp;
  t_sample *imag1 = denominator_spectrum.imagp;
  t_sample *real2 = filter_spectrum.realp;
  t_sample *imag2 = filter_spectrum.imagp;

  t_float filter_amp;

  t_uint from = format == SPECTRUM_REAL ? 1 : 0;
  t_uint to = format == SPECTRUM_REAL ? fft_size >> 1 : (fft_size >> 1) + 1;
  t_uint i;

  if(format == SPECTRUM_REAL) {
    filter_amp = 1. / (1. + beta_in[0] / (real1[0] * real1[0]));
    real2[0] = filter_amp * filter_amp;

    filter_amp = 1. / (1. + beta_in[fft_size >> 1] / (imag1[0] * imag1[0]));
    real2[fft_size >> 1] = filter_amp * filter_amp;
  }

  for(i=from; i<to; i++) {
    filter_amp = 1. / (1. + beta_in[i] / (real1[i] * real1[i] +
					  imag1[i] * imag1[i]));
    real2[i] = filter_amp * filter_amp;
  }

  for(i=(fft_size >> 1)+1; i<fft_size; i++)
    real2[i] = real2[fft_size-i];

  for(i=0; i<fft_size; i++)
    imag2[i] = 0.;

  variable_phase_from_power_spectrum(fft_setup, filter_spectrum, fft_size,
				     phase, true);

  if(format == SPECTRUM_REAL)
    imag2[0] = real2[fft_size >> 1];
}

void make_clip_filter(FFT_Setup *fft_setup, FFT_Split denominator_spectrum, FFT_Split filter_spectrum, t_sample *clip_min, t_sample *clip_max, t_uint fft_size, t_spectrum_format format, double phase)
{
  t_sample *real1 = denominator_spectrum.realp;
  t_sample *imag1 = denominator_spectrum.imagp;
  t_sample *real2 = filter_spectrum.realp;
  t_sample *imag2 = filter_spectrum.imagp;

  t_uint from = format == SPECTRUM_REAL ? 1 : 0;
  t_uint to = format == SPECTRUM_REAL ? fft_size >> 1 : (fft_size >> 1) + 1;
  t_uint i;

  double divisor_power;
  double filter_power;

  if(format == SPECTRUM_REAL) {
    divisor_power = (real1[0] * real1[0]);
    filter_power = divisor_power < clip_min[0] ? divisor_power / clip_min[0] :
      1.;
    filter_power = divisor_power > clip_max[0] ? divisor_power / clip_max[0] :
      filter_power;
    real2[0] = filter_power;

    divisor_power = (imag1[0] * imag1[0]);
    filter_power = divisor_power < clip_min[fft_size >> 1] ? divisor_power /
      clip_min[fft_size >> 1] : 1.;
    filter_power = divisor_power > clip_max[fft_size >> 1] ? divisor_power /
      clip_max[fft_size >> 1] : filter_power;
    real2[fft_size >> 1] = filter_power;
  }

  for(i=from;i<to;i++) {
    divisor_power = (real1[i] * real1[i] + imag1[i] *imag1[i]);
    filter_power = divisor_power < clip_min[i] ? divisor_power / clip_min[i] :
      1.;
    filter_power = divisor_power > clip_max[i] ? divisor_power / clip_max[i] :
      filter_power;
    real2[i] = filter_power;
  }

  for(i=(fft_size >> 1)+1;i<fft_size;i++)
    real2[i] = real2[fft_size-i];

  for(i=0;i<fft_size;i++)
    imag2[i] = 0.;

  variable_phase_from_power_spectrum(fft_setup, filter_spectrum, fft_size,
				     phase, true);

  if(format == SPECTRUM_REAL)
    imag2[0] = real2[fft_size >> 1];
}

void time_to_spectrum(FFT_Setup *fft_setup, t_sample *in_buf, t_uint in_length, FFT_Split spectrum, t_uint fft_size)
{
  t_uint inexact = 0;
  t_uint fft_size_log2 = int_log2(fft_size, &inexact);
  t_uint i;

  if(inexact)
    return;

  for(i=0; i<in_length; i++) {
    spectrum.realp[i] = in_buf[i];
    spectrum.imagp[i] = 0.;
  }

  for(; i<fft_size; i++) {
    spectrum.realp[i] = 0.;
    spectrum.imagp[i] = 0.;
  }

  do_fft(&spectrum, fft_setup, fft_size_log2);   
}

void make_deconvolution_filter(FFT_Setup *fft_setup, FFT_Split denominator_spectrum, FFT_Split filter_spectrum, t_float *filter_specifier, t_float *range_specifier, t_sample filter_db_offset, t_float *filter_in, t_uint filter_length, t_uint fft_size, t_spectrum_format format, t_filter_type mode, double phase, double sample_rate)
{
  t_uint i;

  if(filter_in) {
    time_to_spectrum(fft_setup, filter_in, filter_length,
			   filter_spectrum, fft_size);
    
    if(mode != FILTER_FILTER)
      power_spectrum(filter_spectrum, fft_size, SPECTRUM_FULL);
  }
  else {
    make_freq_dependent_power_array(filter_spectrum.realp, filter_specifier,
    				    fft_size, sample_rate, filter_db_offset);
  }

  switch(mode)
    {
    case FILTER_REGULARISATION:
      make_regularisation_filter(fft_setup, denominator_spectrum,
      				 filter_spectrum, filter_spectrum.realp,
      				 fft_size, format, phase);
      break;
      
    case FILTER_CLIP:
      make_freq_dependent_power_array(filter_spectrum.imagp, range_specifier,
				      fft_size, sample_rate, 0.);
      for(i=0; i<(fft_size>>1); i++)
	filter_spectrum.imagp[i] *= filter_spectrum.realp[i];
      make_clip_filter(fft_setup, denominator_spectrum, filter_spectrum,
		       filter_spectrum.realp, filter_spectrum.imagp, fft_size,
		       format, phase);
      break;

    case FILTER_FILTER:
      if(!filter_in) {
	for(i=0; i<fft_size; i++)
	  filter_spectrum.imagp[i] = 0.;
	variable_phase_from_power_spectrum(fft_setup, filter_spectrum,
					   fft_size, phase, true);
      }

      if(format == SPECTRUM_REAL)
	filter_spectrum.imagp[0] = filter_spectrum.realp[fft_size >> 1];
      break;
    }
}

void deconvolve_variable_phase(FFT_Setup *fft_setup, FFT_Split spectrum_1, FFT_Split spectrum_2, FFT_Split filter_spectrum, t_sample *filter_specifier, t_sample *range_specifier, t_float filter_db_offset, t_float *filter_in, t_uint filter_length, t_uint fft_size, t_spectrum_format format, t_filter_type mode, t_float phase, t_float sample_rate)
{
  make_deconvolution_filter(fft_setup, spectrum_2, filter_spectrum, filter_specifier, range_specifier, filter_db_offset, filter_in, filter_length, fft_size, format, mode, phase, sample_rate);
  deconvolve_with_filter(spectrum_1, spectrum_2, filter_spectrum, fft_size, format);
}

void deconvolve(FFT_Setup *fft_setup, FFT_Split spectrum_1, FFT_Split spectrum_2, FFT_Split filter_spectrum, t_float *filter_specifier, t_float *range_specifier, t_sample filter_db_offset, t_float *filter_in, t_uint filter_length, t_uint fft_size, t_spectrum_format format, t_filter_type mode, t_float phase, t_float delay, t_float sample_rate)
{
  if(phase == 0. && !filter_in) {
    deconvolve_zero_phase(spectrum_1, spectrum_2, filter_spectrum, filter_specifier, range_specifier, filter_db_offset, fft_size, format, mode, sample_rate);
  }
  else {
    deconvolve_variable_phase(fft_setup, spectrum_1, spectrum_2, filter_spectrum, filter_specifier, range_specifier, filter_db_offset, filter_in, filter_length, fft_size, format, mode, phase, sample_rate);
  }

  if(delay)
    delay_spectrum(spectrum_1, fft_size, format, delay);
}


void time_to_halfspectrum(FFT_Setup *fft_setup, t_sample *in_buf, t_uint in_length, FFT_Split spectrum, t_uint fft_size)
{
  t_uint i;
  t_uint inexact = 0;
  t_uint fft_size_log2 = int_log2(fft_size, &inexact);

  if(inexact)
    return;

  unzip_zero(in_buf, &spectrum, in_length, fft_size_log2);
  do_real_fft(&spectrum, fft_setup, fft_size_log2);

  for(i=0; i<(fft_size>>1); i++) {
    spectrum.realp[i] *= 0.5;
    spectrum.imagp[i] *= 0.5;
  }
}

void convolve(FFT_Split fft_data_1, FFT_Split fft_data_2, t_uint fft_size,
	      t_spectrum_format format)
{
  t_sample *real1 = fft_data_1.realp;
  t_sample *real2 = fft_data_2.realp;
  t_sample *imag1 = fft_data_1.imagp;
  t_sample *imag2 = fft_data_2.imagp;

  t_sample a, b, c, d;

  t_uint from = format == SPECTRUM_REAL ? 1 : 0;
  t_uint to = format == SPECTRUM_REAL ? fft_size >> 1 : fft_size;
  t_uint i;

  if(format == SPECTRUM_REAL) {
    real1[0] = real1[0] * real2[0];
    imag1[0] = imag1[0] * imag2[0];
  }

  for(i=from; i<to; i++) {
    a = real1[i];
    b = imag1[i];
    c = real2[i];
    d = imag2[i];

    real1[i] = (a * c) - (b * d);
    imag1[i] = (b * c) + (a * d);
  }
}

void power_full_spectrum_from_half_spectrum(FFT_Split spectrum, t_uint fft_size)
{
  t_uint i;
  
  t_sample *real = spectrum.realp;
  t_sample *imag = spectrum.imagp;

  real[0] *= spectrum.realp[0];

  for(i=1; i<fft_size>>1; i++)
    real[i] = (real[i] * real[i]) + (imag[i] * imag[i]);

  real[i++] = imag[0] * imag[0];

  for(; i<fft_size; i++)
    real[i] = real[fft_size-i];

  for(i=0; i<fft_size; i++)
    imag[i] = 0.;
}

void spectrum_to_time(FFT_Setup *fft_setup, t_sample *out_buf, FFT_Split spectrum, t_uint fft_size, t_spectrum_format half_spectrum)
{
  t_uint inexact = 0;
  t_uint fft_size_log2 = int_log2(fft_size, &inexact);
  t_uint i = 0;

  double scale = 1. / (double)fft_size;

  if(inexact)
    return;

  if(half_spectrum == SPECTRUM_FULL) {
    spectrum.imagp[0] = spectrum.realp[fft_size >> 1];
  }

  do_real_ifft(&spectrum, fft_setup, fft_size_log2);
  zip_sample(&spectrum, out_buf, (t_uint) 1 << (fft_size_log2 - (t_uint) 1));

  for(i=0; i<fft_size; i++)
    out_buf[i] *= scale;
}



void setup_hann_wind()
{
  t_uint i;
  for(i=0; i<4097; i++)
    hann_table[i] = 0.5 * (1. - cos(M_PI * (4095 - i) / (t_float) 4095));
}

t_float fast_hann_wind(t_float in)
{
  t_uint index;

  t_float lo, hi, f_index, fract;

  f_index = in * 4095;
  index = (t_uint)f_index;
  fract = f_index - index;

  lo = hann_table[index];
  hi = hann_table[index+1];

  return lo + fract * (hi - lo);
}

void smooth_power_spectrum(FFT_Split spectrum, long mode, t_uint fft_size, t_float smooth_lo, t_float smooth_hi)
{
  t_sample *spectrum_out = spectrum.realp;
  t_sample *spectrum_in = spectrum.imagp;

  t_float filter, filter_val, half_width_recip, oct_width, accum, smooth_mul;

  t_int lo, hi, left, half_width;
  t_int nyquist_bin = (fft_size >> 1);
  t_int fft_size_m1 = fft_size - 1;
  t_int limit = fft_size_m1;
  t_int i, j;

  /*  probably there is a danger with the upper branch if the smooth hi is 1 */
  
  smooth_lo = smooth_lo > 1. ? 1. : smooth_lo;
  smooth_hi = smooth_hi > 1. ? 1. : smooth_hi;
  smooth_mul = smooth_hi - smooth_lo;

  for(i=0; i<(t_int)fft_size; i++)
    spectrum_in[i] = spectrum_out[i];

  switch(mode)
    {
    case 1:

      if(!hann_setup_flag) {
	setup_hann_wind();
	hann_setup_flag = 1;
      }

      for(i=0; i<nyquist_bin+1; i++) {
	half_width = (t_int) (((( (t_float)i / (t_float)nyquist_bin) * smooth_mul) + smooth_lo) * (t_float)nyquist_bin);
	filter_val = spectrum_in[i];
	half_width_recip = half_width > 1. ? 2. / (2 * half_width - 1.) : 1.;

	for(j=1; j<half_width; j++) {
	  left = (i - j) & fft_size_m1;
	  filter = fast_hann_wind(j * half_width_recip);
	  filter_val += filter * (spectrum_in[left] + spectrum_in[i + j]);
	}

	spectrum_out[i] = filter_val * half_width_recip;
      }
      break;

    case 2:

      for(i=1; i<(t_int)fft_size; i++)
	spectrum_in[i] += spectrum_in[i - 1];

      for(i=0; i<nyquist_bin+1; i++) {
	half_width = (t_int) (((( (t_float)i / (t_float)nyquist_bin) * smooth_mul) + smooth_lo) * (t_float) nyquist_bin);
	accum = 0.;

	lo = i - half_width;
	hi = i + half_width;

	if(lo<0) {
	  accum = (spectrum_in[-lo] - spectrum_in[0]) / -lo;
	  lo = 0;
	}

	if(lo == hi)
	  hi++;

	if(hi>limit)
	  hi = limit;

	spectrum_out[i] = accum + (spectrum_in[hi] - spectrum_in[lo]) / (hi - lo);
	
      }
      
      break;
    case 3:

      for(i=1; i<(t_int)fft_size; i++)
	spectrum_in[i] += spectrum_in[i - 1];

      spectrum_out[0] = spectrum_in[0];

      for(i=1; i<nyquist_bin+1; i++) {
	oct_width = ((( (t_float)i / (t_float)nyquist_bin) * smooth_mul) + smooth_lo);
	oct_width = pow(2., oct_width * 0.5);

	lo = (t_int)(i/oct_width);
	hi = (t_int)(i*oct_width);

	if(lo==hi)
	  lo--;

	if(hi>limit)
	  hi = limit;

	spectrum_out[i] = (spectrum_in[hi] - spectrum_in[lo]) / (hi-lo);
      }
      
      break;
    }

  for(i=nyquist_bin+1; i<(t_int)fft_size; i++)
    spectrum_out[i] = spectrum_out[fft_size - i];

  for(i=0; i<(t_int)fft_size; i++)
    spectrum_in[i] = 0.;
}
