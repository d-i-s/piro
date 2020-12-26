#include "m_pd.h"

#include "./lib/u_multi_channel_convolve.h"
#include "./lib/z_array.h"

#define MAXIMUM_DSP_CHANS 64

static t_class *multiconvolve_tilde_class;

typedef struct _multiconvolve_tilde {
  t_object  x_obj;
  t_sample x_f;

  t_multi_channel_convolve *multi;

  long parallel_mode;
  long num_in_chans;
  long num_out_chans;

  t_sample *ins[MAXIMUM_DSP_CHANS];
  t_sample *outs[MAXIMUM_DSP_CHANS];
	
  long fixed_impulse_length;
  
} t_multiconvolve_tilde;


static void multiconvolve_tilde_fixed_size_set(t_multiconvolve_tilde *x, t_symbol *sym, long argc, t_atom *argv)
{
  long i, j;
  long error = 0;

  if(argc) {

    long new_fixed_size = (long)atom_getint(argv);

    new_fixed_size = new_fixed_size < 0 ? 0 : new_fixed_size;
		
    if (new_fixed_size != x->fixed_impulse_length && x->multi) {
      if (!x->parallel_mode) {
        for (i = 0; i < x->num_out_chans; i++)
	  for (j = 0; j < x->num_in_chans; j++)
	    if (multi_channel_convolve_resize(x->multi, j, i, new_fixed_size ? (t_uint) new_fixed_size : 16384UL) == CONVOLVE_ERR_MEM_UNAVAILABLE)
	      error = 1;
      }
      else {
        for (i = 0; i < x->num_in_chans; i++)
	  if (multi_channel_convolve_resize(x->multi, i, i, new_fixed_size ? (t_uint) new_fixed_size : 16384UL) == CONVOLVE_ERR_MEM_UNAVAILABLE)
	    error = 1;
      }
			
      if (error)
        pd_error(x, "could not allocate memory for fixed size");
    }
    x->fixed_impulse_length = (long)new_fixed_size;
  }
	
}

static void multiconvolve_tilde_clear(t_multiconvolve_tilde *x)
{	
  if (!x->multi)
    return;
	
  multi_channel_convolve_clear(x->multi, x->fixed_impulse_length ? false : true);
}

static void multiconvolve_tilde_set(t_multiconvolve_tilde *x, t_symbol *sym, long argc, t_atom *argv)
{
  t_convolve_error set_error = CONVOLVE_ERR_NONE;

  t_symbol *buffer;
  t_float *temp;
	
  long in_chan;
  long out_chan;
	
  t_int impulse_length;
	
  if (argc < 1)	{
    pd_error(x, "not enough arguments to message %s", sym->s_name);
    return;
  }
	
  /* Get input channel if present */
	
  if (argv->a_type == A_FLOAT) {
    in_chan = (long)atom_getint(argv++) - 1;
    argc--;
  }		
  else
    in_chan = 0;

  /* Get out channel if present */
	
  if (argv->a_type == A_FLOAT) {
    if (argc == 1) {
      pd_error(x, "no buffer given for set message");
      return;
    }
    out_chan = (long)atom_getint(argv++) - 1;
    argc--;
  }
  else
    out_chan = in_chan;
	
  if (x->parallel_mode && out_chan != in_chan) {
    pd_error((t_object *)x, "parallel mode disallows mismatching input/output channels (%ld to %ld requested)", in_chan + 1, out_chan + 1);
    return;
  }
	
  /* Get buffer */
	
  buffer = atom_getsymbol(argv++);
  argc--;
	
  /* Check buffer */
	
  if(!buffer_check(buffer)) {
    pd_error(x, "array not found %s", buffer->s_name);
    return;
  }
	
  /* Load to temporary buffer */
	
  impulse_length = buffer_length(buffer);
  temp = (t_float *)aligned_getbytes(impulse_length * sizeof(t_float));
  
  if (!temp) {
    pd_error(x, "could not allocate temporary memory for temporary storage");
    return;
  }
		
  buffer_read(buffer, temp, impulse_length);
  
  if(x->multi)
    set_error = multi_channel_convolve_set(x->multi, (t_uint) in_chan, (t_uint) out_chan, temp, (t_uint)impulse_length, x->fixed_impulse_length ? false : true);
	
  if (set_error == CONVOLVE_ERR_IN_CHAN_OUT_OF_RANGE)
    pd_error(x, "input channel %ld out of range", in_chan + 1);
  if (set_error == CONVOLVE_ERR_OUT_CHAN_OUT_OF_RANGE)
    pd_error(x, "output channel %ld out of range", out_chan + 1);
  if (set_error == CONVOLVE_ERR_MEM_UNAVAILABLE)		
    pd_error(x, "memory unavailable for set operation");
  if (set_error == CONVOLVE_ERR_MEM_ALLOC_TOO_SMALL)		
    pd_error(x, "requested buffer / length too large for set fixed size");

  
  aligned_freebytes(temp, impulse_length * sizeof(t_float));
}


static t_int *multiconvolve_tilde_perform(t_int *w)
{
  t_multiconvolve_tilde *x = (t_multiconvolve_tilde *)(w[1]);
  long vec_size = (long)(w[2]);


  multi_channel_convolve_process(x->multi, x->ins, x->outs, 0, 0, (t_uint)vec_size, (t_uint)x->num_in_chans, (t_uint)x->num_out_chans);

  
  return (w+3);
}

static void multiconvolve_tilde_dsp(t_multiconvolve_tilde *x, t_signal **sp)
{
  long i;
  

  for(i=0; i<x->num_in_chans; i++)
    x->ins[i] = sp[i]->s_vec;

  for(i=0; i<x->num_out_chans; i++)
    x->outs[i] = sp[i+x->num_in_chans]->s_vec;


  if(x->multi)
    dsp_add(multiconvolve_tilde_perform, 2, x, sp[0]->s_n);
}

static void multiconvolve_tilde_free(t_multiconvolve_tilde *x)
{
  multi_channel_convolve_free(x->multi);
}

static void *multiconvolve_tilde_new(t_symbol *s, short argc, t_atom *argv)
{
  t_multiconvolve_tilde *x = (t_multiconvolve_tilde *)pd_new(multiconvolve_tilde_class);

  t_int num_out_chans = 1;
  t_int num_in_chans = 1;
  t_int latency_mode = 0;

  /* long actual_num_out_chans; */
  long i;

  if(!x)
    return NULL;

  if(argc && argv->a_type != A_SYMBOL) {
    num_in_chans = atom_getint(argv++);
    num_in_chans = num_in_chans < 1 ? 1 : num_in_chans;
    num_in_chans = num_in_chans > MAXIMUM_DSP_CHANS ?
      MAXIMUM_DSP_CHANS : num_in_chans;
    argc--;
  }

  if(argc && argv->a_type != A_SYMBOL) {
    num_out_chans = atom_getint(argv++);
    num_out_chans = num_out_chans < 1 ? 1 : num_out_chans;
    num_out_chans = num_out_chans > MAXIMUM_DSP_CHANS ?
      MAXIMUM_DSP_CHANS : num_out_chans;
  }

  if(argc && argv->a_type == A_SYMBOL && atom_getsymbol(argv) != gensym("fixedsize")) {
    t_symbol *mode = atom_getsymbol(argv++);

    if(mode == gensym("short"))
      latency_mode = 1;
    if(mode == gensym("medium"))
      latency_mode = 2;

    if(!latency_mode && mode != gensym("zero"))
      pd_error(x, "unknown latency mode - %s", mode->s_name);

    argc--;
  }

  /* actual_num_out_chans = (long) (num_out_chans <= 0 ? num_in_chans : num_out_chans); */
 
  x->num_in_chans = (long)num_in_chans;
  x->num_out_chans = (long)num_out_chans;
  x->fixed_impulse_length = 0;

  if(!num_out_chans)
    x->parallel_mode = 1;
  else
    x->parallel_mode = 0;

  for(i=1;i<x->num_in_chans;i++)
    inlet_new(&x->x_obj,&x->x_obj.ob_pd,&s_signal,&s_signal);
  
  for(i=0;i<x->num_out_chans;i++)
    outlet_new(&x->x_obj, &s_signal);

  x->multi = multi_channel_convolve_new((t_uint)num_in_chans, (t_uint)num_out_chans, (t_convolve_latency_mode)latency_mode, (t_uint)16384);

  return (void *)x;
}

void multiconvolve_tilde_setup(void)
{
#if PD_FLOATSIZE == 32
  multiconvolve_tilde_class = class_new(gensym("multiconvolve~"),
					(t_newmethod)multiconvolve_tilde_new,
					(t_method)multiconvolve_tilde_free,
					sizeof(t_multiconvolve_tilde),
					CLASS_DEFAULT,
					A_GIMME, 0);
#elif PD_FLOATSIZE == 64
  multiconvolve_tilde_class = class_new64(gensym("multiconvolve~"),
					  (t_newmethod)multiconvolve_tilde_new,
					  (t_method)multiconvolve_tilde_free,
					  sizeof(t_multiconvolve_tilde),
					  CLASS_DEFAULT, 
					  A_GIMME, 0);
#else
#error [multiconvolve~]: invalid FLOATSIZE: must be 32 or 64
#endif
  CLASS_MAINSIGNALIN(multiconvolve_tilde_class, t_multiconvolve_tilde, x_f);
  class_addmethod(multiconvolve_tilde_class,
		  (t_method)multiconvolve_tilde_dsp, gensym("dsp"), 0);
  
  class_addmethod(multiconvolve_tilde_class,
		  (t_method)multiconvolve_tilde_set, gensym("set"), A_GIMME, 0);
  class_addmethod(multiconvolve_tilde_class,
		  (t_method)multiconvolve_tilde_clear, gensym("clear"), 0);
  class_addmethod(multiconvolve_tilde_class,
		  (t_method)multiconvolve_tilde_fixed_size_set,
		  gensym("fixedsize"), A_GIMME, 0);
}
