
/*
 *  routine to access the array
 *	
 *
 *  Copyright 2019/20 Marco Matteo Markidis. All rights reserved.
 *
 */


#include "z_array.h"

int buffer_check(t_symbol *name)
{
  int found = 0;

  t_garray *a = NULL;

  if(!(a = (t_garray *)pd_findbyclass(name, garray_class)))
    return found;
  else
    found = 1;

  return found;
}

int attach_array(t_symbol *x_name, t_garray **x_buf, t_word **x_samples,
		 int *x_frames)
{
  t_garray *a = NULL;
  t_symbol *name = x_name;
  t_word *samples;
  int frames;
  int err = 0;

  if(!name) {
    return err;
  }
  
  if(!(a = (t_garray *)pd_findbyclass(name, garray_class))) {
    x_samples = 0;
    if(name->s_name) {
      post("%s: no such array", name->s_name);
      return err;
    }
  }

  if(!garray_getfloatwords(a, &frames, &samples)) {
    x_samples = 0;
    post("bad array for %s", name->s_name);
    return err;
  }
  else {
    *x_frames = frames;
    *x_samples = samples;
    *x_buf = a;
    return 1;
  }
}

int buffer_multiple_names(t_symbol *sin[], t_symbol *sout[], t_int length[], int argc, t_atom *argv, t_int in_place, t_int *overall_length, t_int *max_length)
{
  int err = 0;
  int num_buffers = 0;

  t_garray *buf;
  t_symbol *buf_name;
  t_word *buf_samples;
  int buf_frames = 0;

  t_int local_max = 0;
  t_int local_overall = 0;

  int i;

  num_buffers = in_place ? argc : argc >> 1;

  if(in_place) {
    for(i=0; i<num_buffers; i++) {
      buf_name = atom_getsymbol(argv++);
      err = attach_array(buf_name, &buf, &buf_samples, &buf_frames);
      if(!err)
	return err;
      else {
	if(sin)
	  sin[i] = buf_name;
	if(sout)
	  sout[i] = buf_name;
	if(length)
	  length[i] = buf_frames;
	local_overall += buf_frames;
	local_max = buf_frames > local_max ? buf_frames : local_max;
      }
    }
  }
  else {
    for(i=0; i<num_buffers; i++) {
      buf_name = atom_getsymbol(argv++);
      err = attach_array(buf_name, &buf, &buf_samples, &buf_frames);
      if(!err)
	return err;
      else {
	sin[i] = buf_name;
	length[i] = buf_frames;
	local_overall += buf_frames;
	local_max = buf_frames > local_max ? buf_frames : local_max;
      }
    }
    for(; i<argc; i++) {
      buf_name = atom_getsymbol(argv++);
      err = attach_array(buf_name, &buf, &buf_samples, &buf_frames);
      if(!err)
	return err;
      else {
	sout[i-num_buffers] = buf_name;
	length[i] = buf_frames;
	local_overall += buf_frames;
	local_max = buf_frames > local_max ? buf_frames : local_max;
      }
    }
  }

  *overall_length = local_overall;
  *max_length = local_max;
  
  return num_buffers;
}

int buffer_read(t_symbol *s, t_float *out, t_int total_length)
{
  int err = 0;
  t_int i;

  t_garray *buf;
  t_symbol *buf_name;
  t_word *buf_samples = NULL;
  int buf_frames = 0;

  buf_name = s;
  
  if(buf_name)
    err = attach_array(buf_name, &buf, &buf_samples, &buf_frames);
  else
    return err;
  
  if(!err)
    return err;

  for(i=0; i<buf_frames; i++)
    *out++ = buf_samples[i].w_float;

  for(; i<total_length; i++)
    *out++ = 0.;

  return buf_frames;
}

int buffer_length(t_symbol *s)
{
  t_garray *buf;
  t_symbol *buf_name;
  t_word *buf_samples;
  int buf_frames = 0;

  int err = 0;

  buf_name = s;

  if(buf_name)
    err = attach_array(buf_name, &buf, &buf_samples, &buf_frames);
  else
    return err;
  
  if(!err)
    return err;

  return buf_frames;
}

int buffer_write(t_symbol *s, t_sample *in, t_int write_length, t_int resize,
		 t_float mul)
{
  t_int i;

  t_garray *buf;
  t_symbol *buf_name;
  t_word *buf_samples;
  int buf_frames = 0;

  int err = 0;

  buf_name = s;
  
  err = attach_array(buf_name, &buf, &buf_samples, &buf_frames);
  
  if(!err)
    {
      post("error in writing buffer");
      return err;
    }
  
  if(resize) {
    garray_resize_long(buf, write_length);
    err = attach_array(buf_name, &buf, &buf_samples, &buf_frames); 
  if(!err)
    {
      post("error in writing buffer");
      return err;
    }
  }
  
  for(i=0; i<write_length; i++)
    {
      t_sample f = (t_float)(in[i] * mul);//*in++ * mul;
      if(PD_BIGORSMALL(f))
	f = 0;
      buf_samples[i].w_float = f;
      //(buf_samples)->w_float = f; // attento, ++
    
    }
  //buf_samples[i].w_float = (t_float)(in[i] * mul);
  // (buf_samples++)->w_float = (t_float)(*in++ * mul);

  if(!buf_samples) post("error in writing buffer");

  //buf_samples = vec;
  
  if(resize==0) {
    for(;i<buf_frames;i++) {
      buf_samples[i].w_float = 0.;
      //(buf_samples++)->w_float = 0.;
    }
  }

  garray_redraw(buf);
  
  return 1;
}
