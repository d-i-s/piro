#include "z_utils.h"

t_float delay_retriever(t_atom *delay, t_uint fft_size, t_float sample_rate)
{
  return (delay->a_type == A_SYMBOL) ? (fft_size >> 1) : (atom_getfloat(delay)*sample_rate/1000.);
}

t_symbol *filter_retriever(t_atom *specifier)
{
  if(specifier->a_type == A_SYMBOL)
    return specifier->a_w.w_symbol;
  else 
    return NULL;
}

t_float phase_retriever(t_atom a)
{
  t_symbol *sym = atom_getsymbol(&a);

  if(a.a_type == A_FLOAT)
    return a.a_type;

  if(sym == gensym("lin") || sym == gensym("linear"))
    return 0.5;

  if(sym==gensym("min") || sym == gensym("minimum"))
    return 0.;

  if(sym==gensym("max") || sym == gensym("maximum"))
    return 1.;

  return 0;
}

t_atom phase_parser(t_atom a)
{
  if(a.a_type == A_FLOAT) {
    double phase = atom_getfloat(&a);
    phase = phase < 0. ? 0. : phase;
    phase = phase > 1. ? 1. : phase;

    SETFLOAT(&a, phase);
  }
  else {
    t_symbol *sym = atom_getsymbol(&a);

    if(sym != gensym("lin") && sym != gensym("linear") && sym != gensym("min") && sym != gensym("minimum") && sym != gensym("max") && sym != gensym("maximum"))
      SETSYMBOL(&a, gensym("minimum"));
  }

  return a;
}
