#ifndef __Z_ARRAY__
#define __Z_ARRAY__

#include "m_pd.h"

int buffer_check(t_symbol *name);
int attach_array(t_symbol *x_name, t_garray **x_buf, t_word **x_samples,
		 int *x_frames);
int buffer_multiple_names(t_symbol *sin[], t_symbol *sout[], t_int length[], int argc, t_atom *argv, t_int in_place, t_int *overall_length, t_int *max_length);
int buffer_read(t_symbol *s, t_float *out, t_int total_length);
int buffer_length(t_symbol *s);
int buffer_write(t_symbol *s, t_sample *in, t_int write_length, t_int resize,
		 t_float mul);


#endif /* __Z_ARRAY__ */
