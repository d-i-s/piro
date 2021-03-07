#ifndef __ZUTILS__
#define __ZUTILS__

//#include "m_pd.h"
#include "z_fft.h"

t_float delay_retriever(t_atom *delay, t_uint fft_size, t_float sample_rate);
t_symbol *filter_retriever(t_atom *specifier);
t_float phase_retriever(t_atom a);
t_atom phase_parser(t_atom a);

#endif /* __ZUTILS__ */
