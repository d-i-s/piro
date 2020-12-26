#include "z_fft.h"

static int SSE2_check()
{
#if defined (UNIX) || defined ( MACOSX ) 
  return 1;
#elif defined (NT)
  int SSE2_flag = 0;
  int CPUInfo[4] = {-1, 0, 0, 0};
  int nIds;

  __cpuid(CPUInfo, 0);
  nIds = CPUInfo[0];

  if(nIds > 0) {
    __cpuid(CPUInfo, 1);
    SSE2_flag = (CPUInfo[3] >> 26) & 0x1;
  }

  return SSE2_flag;
#endif
}

void pass_1_2_reorder_simd32(FFT_Split32 *input, t_uint length)
{
  t_uint i;
	
  vFloat *r1_ptr = (vFloat *)input->realp;
  vFloat *r2_ptr = r1_ptr + (length >> 4);
  vFloat *r3_ptr = r2_ptr + (length >> 4);
  vFloat *r4_ptr = r3_ptr + (length >> 4);
	
  vFloat *i1_ptr = (vFloat *)input->imagp;
  vFloat *i2_ptr = i1_ptr + (length >> 4);
  vFloat *i3_ptr = i2_ptr + (length >> 4);
  vFloat *i4_ptr = i3_ptr + (length >> 4);
	
  vFloat r1, r2, r3, r4, r5, r6, r7, r8;
  vFloat i1, i2, i3, i4, i5, i6, i7, i8;
	
  for (i = 0; i < length >> 4; i++)
    {
      r5 = *r1_ptr;
      r6 = *r2_ptr;
      r3 = *r3_ptr;
      r4 = *r4_ptr;
		
      i5 = *i1_ptr;
      i6 = *i2_ptr;
      i3 = *i3_ptr;
      i4 = *i4_ptr;
		
      r1 = F32_VEC_ADD_OP(r5, r3);
      r2 = F32_VEC_ADD_OP(r6, r4);
      r3 = F32_VEC_SUB_OP(r5, r3);
      r4 = F32_VEC_SUB_OP(r6, r4);
		
      i1 = F32_VEC_ADD_OP(i5, i3);
      i2 = F32_VEC_ADD_OP(i6, i4);
      i3 = F32_VEC_SUB_OP(i5, i3);
      i4 = F32_VEC_SUB_OP(i6, i4);
		
      r5 = F32_VEC_ADD_OP(r1, r2);
      r6 = F32_VEC_SUB_OP(r1, r2);
      r7 = F32_VEC_ADD_OP(r3, i4);
      r8 = F32_VEC_SUB_OP(r3, i4);
		
      i5 = F32_VEC_ADD_OP(i1, i2);
      i6 = F32_VEC_SUB_OP(i1, i2);
      i7 = F32_VEC_SUB_OP(i3, r4);
      i8 = F32_VEC_ADD_OP(i3, r4);
		
      r1 = F32_VEC_SHUFFLE(r5, r7, F32_SHUFFLE_CONST(1, 0, 1, 0));
      r2 = F32_VEC_SHUFFLE(r5, r7, F32_SHUFFLE_CONST(3, 2, 3, 2));
      r3 = F32_VEC_SHUFFLE(r6, r8, F32_SHUFFLE_CONST(1, 0, 1, 0));
      r4 = F32_VEC_SHUFFLE(r6, r8, F32_SHUFFLE_CONST(3, 2, 3, 2));
		
      *r1_ptr++ = F32_VEC_SHUFFLE(r1, r3, F32_SHUFFLE_CONST(2, 0, 2, 0));
      *r2_ptr++ = F32_VEC_SHUFFLE(r2, r4, F32_SHUFFLE_CONST(2, 0, 2, 0));
      *r3_ptr++ = F32_VEC_SHUFFLE(r1, r3, F32_SHUFFLE_CONST(3, 1, 3, 1));
      *r4_ptr++ = F32_VEC_SHUFFLE(r2, r4, F32_SHUFFLE_CONST(3, 1, 3, 1));
		
      i1 = F32_VEC_SHUFFLE(i5, i7, F32_SHUFFLE_CONST(1, 0, 1, 0));
      i2 = F32_VEC_SHUFFLE(i5, i7, F32_SHUFFLE_CONST(3, 2, 3, 2));
      i3 = F32_VEC_SHUFFLE(i6, i8, F32_SHUFFLE_CONST(1, 0, 1, 0));
      i4 = F32_VEC_SHUFFLE(i6, i8, F32_SHUFFLE_CONST(3, 2, 3, 2));
		
      *i1_ptr++ = F32_VEC_SHUFFLE(i1, i3, F32_SHUFFLE_CONST(2, 0, 2, 0));
      *i2_ptr++ = F32_VEC_SHUFFLE(i2, i4, F32_SHUFFLE_CONST(2, 0, 2, 0));
      *i3_ptr++ = F32_VEC_SHUFFLE(i1, i3, F32_SHUFFLE_CONST(3, 1, 3, 1));
      *i4_ptr++ = F32_VEC_SHUFFLE(i2, i4, F32_SHUFFLE_CONST(3, 1, 3, 1));
    }
}

void pass_3_reorder_simd32(FFT_Split32 *input, t_uint length)
{
  t_uint offset = length >> 5;
  t_uint outer_loop = length >> 6;
  t_uint i, j;
	
  vFloat r1, r2, r3, r4, r5;
  vFloat i1, i2, i3, i4, i5;
  vFloat twiddle_c = {1, (t_sample) SQRT_2_2, 0, (t_sample) -SQRT_2_2};
  vFloat twiddle_s = {0, (t_sample)-SQRT_2_2, -1, (t_sample) -SQRT_2_2};
	
  vFloat *r1_ptr = (vFloat *)input->realp;
  vFloat *i1_ptr = (vFloat *)input->imagp;
  vFloat *r2_ptr = r1_ptr + offset;
  vFloat *i2_ptr = i1_ptr + offset;
	
  for (i = 0, j = 0; i < length >> 1; i += 8)
    {
      // Get input
		
      r1 = *r1_ptr;
      i1 = *i1_ptr;
      r3 = *(r1_ptr + 1);
      i3 = *(i1_ptr + 1);
      r2 = *r2_ptr;
      i2 = *i2_ptr;
      r4 = *(r2_ptr + 1);
      i4 = *(i2_ptr + 1);
		
      // Multiply by twiddle
		
      r5 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r2, twiddle_c),
			  F32_VEC_MUL_OP(i2, twiddle_s));
      i5 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r2, twiddle_s),
			  F32_VEC_MUL_OP(i2, twiddle_c));
		
      // Store output (swapping as necessary)
		
      *r1_ptr = F32_VEC_ADD_OP(r1, r5);
      *i1_ptr = F32_VEC_ADD_OP(i1, i5);
		
      *(r1_ptr + 1) = F32_VEC_SUB_OP(r1, r5);
      *(i1_ptr + 1) = F32_VEC_SUB_OP(i1, i5);
		
      // Multiply by twiddle
		
      r5 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r4, twiddle_c),
			  F32_VEC_MUL_OP(i4, twiddle_s));
      i5 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r4, twiddle_s),
			  F32_VEC_MUL_OP(i4, twiddle_c));
		
      // Store output (swapping as necessary)
		
      *r2_ptr = F32_VEC_ADD_OP(r3, r5);
      *i2_ptr = F32_VEC_ADD_OP(i3, i5);
		
      *(r2_ptr + 1) = F32_VEC_SUB_OP(r3, r5);
      *(i2_ptr + 1) = F32_VEC_SUB_OP(i3, i5);
		
      r1_ptr += 2;
      r2_ptr += 2;
      i1_ptr += 2;
      i2_ptr += 2;
		
      if (!(++j % outer_loop))
	{
	  r1_ptr += offset;
	  r2_ptr += offset;
	  i1_ptr += offset;
	  i2_ptr += offset;
	}
    }
}

void pass_trig_table_reorder_simd32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint length, t_uint pass)
{
  t_uint size = 2 << pass;
  t_uint incr = size >> 3;
  t_uint loop = size;
  t_uint offset = length >> (pass + 3);
  t_uint outer_loop = ((length >> 1) / size) / ((t_uint) 1 << pass);
  t_uint i, j;
	
  vFloat r1, r2, r3, r4, r5;
  vFloat i1, i2, i3, i4, i5;
  vFloat twiddle_c, twiddle_s;
	
  vFloat *r1_ptr = (vFloat *)input->realp;
  vFloat *i1_ptr = (vFloat *)input->imagp;
  vFloat *r2_ptr = r1_ptr + offset;
  vFloat *i2_ptr = i1_ptr + offset;
	
  for (j = 0, i = 0; i < length >> 1; loop += size)
    {
      vFloat *table_r = (vFloat *)setup->tables[pass - PASS_TRIG_OFFSET].realp;
      vFloat *table_i = (vFloat *)setup->tables[pass - PASS_TRIG_OFFSET].imagp;
		
      for (; i < loop; i += 8)
	{
	  // Get input
			
	  r1 = *r1_ptr;
	  i1 = *i1_ptr;
	  r2 = *r2_ptr;
	  i2 = *i2_ptr;
			
	  // Get Twiddle
			
	  twiddle_c = *table_r++;
	  twiddle_s = *table_i++;
			
	  // Multiply by twiddle
			
	  r5 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r2, twiddle_c),
			      F32_VEC_MUL_OP(i2, twiddle_s));
	  i5 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r2, twiddle_s),
			      F32_VEC_MUL_OP(i2, twiddle_c));
			
	  // Get input
			
	  r3 = *(r1_ptr + incr);
	  i3 = *(i1_ptr + incr);
	  r4 = *(r2_ptr + incr);
	  i4 = *(i2_ptr + incr);
			
	  // Store output (swapping as necessary)
			
	  *r1_ptr = F32_VEC_ADD_OP(r1, r5);
	  *i1_ptr = F32_VEC_ADD_OP(i1, i5);
			
	  *(r1_ptr++ + incr) = F32_VEC_SUB_OP(r1, r5);
	  *(i1_ptr++ + incr) = F32_VEC_SUB_OP(i1, i5);
			
	  // Multiply by twiddle
			
	  r5 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r4, twiddle_c),
			      F32_VEC_MUL_OP(i4, twiddle_s));
	  i5 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r4, twiddle_s),
			      F32_VEC_MUL_OP(i4, twiddle_c));
			
	  // Store output (swapping as necessary)
			
	  *r2_ptr = F32_VEC_ADD_OP(r3, r5);
	  *i2_ptr = F32_VEC_ADD_OP(i3, i5);
			
	  *(r2_ptr++ + incr) = F32_VEC_SUB_OP(r3, r5);
	  *(i2_ptr++ + incr) = F32_VEC_SUB_OP(i3, i5);
	}
		
      r1_ptr += incr;
      r2_ptr += incr;
      i1_ptr += incr;
      i2_ptr += incr;
		
      if (!(++j % outer_loop))
	{
	  r1_ptr += offset;
	  r2_ptr += offset;
	  i1_ptr += offset;
	  i2_ptr += offset;
	}
    }
}

void pass_trig_table_simd32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint length, t_uint pass)
{
  
  t_uint size = 2 << pass;
  t_uint incr = size >> 3;
  t_uint loop = size;
  t_uint i;
	
  vFloat r1, r2, r3, i1, i2, i3;
  vFloat twiddle_c, twiddle_s;
	
  vFloat *r1_ptr = (vFloat *)input->realp;
  vFloat *i1_ptr = (vFloat *)input->imagp;
  vFloat *r2_ptr = r1_ptr + incr;
  vFloat *i2_ptr = i1_ptr + incr;
	
  for (i = 0; i < length; loop += size)
    {
      vFloat *table_r = (vFloat *)setup->tables[pass - PASS_TRIG_OFFSET].realp;
      vFloat *table_i = (vFloat *)setup->tables[pass - PASS_TRIG_OFFSET].imagp;
		
      for (; i < loop; i += 8)
	{
	  // Get input
			
	  r1 = *r1_ptr;
	  i1 = *i1_ptr;
	  r2 = *r2_ptr;
	  i2 = *i2_ptr;
			
	  // Get Twiddle
			
	  twiddle_c = *table_r++;
	  twiddle_s = *table_i++;
			
	  // Multiply by twiddle
			
	  r3 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r2, twiddle_c),
			      F32_VEC_MUL_OP(i2, twiddle_s));
	  i3 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r2, twiddle_s),
			      F32_VEC_MUL_OP(i2, twiddle_c));
			
	  // Store output (same pos as inputs)
			
	  *r1_ptr++ = F32_VEC_ADD_OP(r1, r3);
	  *i1_ptr++ = F32_VEC_ADD_OP(i1, i3);
			
	  *r2_ptr++ = F32_VEC_SUB_OP(r1, r3);
	  *i2_ptr++ = F32_VEC_SUB_OP(i1, i3);
	}
		
      r1_ptr += incr;
      r2_ptr += incr;
      i1_ptr += incr;
      i2_ptr += incr;
    }
}

#ifdef PD_SAMPLE_PRECISION /* t_sample == float */
void pass_1_2_reorder_simd(FFT_Split *input, t_uint length)
{
  t_uint i;
	
  vFloat *r1_ptr = (vFloat *)input->realp;
  vFloat *r2_ptr = r1_ptr + (length >> 4);
  vFloat *r3_ptr = r2_ptr + (length >> 4);
  vFloat *r4_ptr = r3_ptr + (length >> 4);
	
  vFloat *i1_ptr = (vFloat *)input->imagp;
  vFloat *i2_ptr = i1_ptr + (length >> 4);
  vFloat *i3_ptr = i2_ptr + (length >> 4);
  vFloat *i4_ptr = i3_ptr + (length >> 4);
	
  vFloat r1, r2, r3, r4, r5, r6, r7, r8;
  vFloat i1, i2, i3, i4, i5, i6, i7, i8;
	
  for (i = 0; i < length >> 4; i++)
    {
      r5 = *r1_ptr;
      r6 = *r2_ptr;
      r3 = *r3_ptr;
      r4 = *r4_ptr;
		
      i5 = *i1_ptr;
      i6 = *i2_ptr;
      i3 = *i3_ptr;
      i4 = *i4_ptr;
		
      r1 = F32_VEC_ADD_OP(r5, r3);
      r2 = F32_VEC_ADD_OP(r6, r4);
      r3 = F32_VEC_SUB_OP(r5, r3);
      r4 = F32_VEC_SUB_OP(r6, r4);
		
      i1 = F32_VEC_ADD_OP(i5, i3);
      i2 = F32_VEC_ADD_OP(i6, i4);
      i3 = F32_VEC_SUB_OP(i5, i3);
      i4 = F32_VEC_SUB_OP(i6, i4);
		
      r5 = F32_VEC_ADD_OP(r1, r2);
      r6 = F32_VEC_SUB_OP(r1, r2);
      r7 = F32_VEC_ADD_OP(r3, i4);
      r8 = F32_VEC_SUB_OP(r3, i4);
		
      i5 = F32_VEC_ADD_OP(i1, i2);
      i6 = F32_VEC_SUB_OP(i1, i2);
      i7 = F32_VEC_SUB_OP(i3, r4);
      i8 = F32_VEC_ADD_OP(i3, r4);
		
      r1 = F32_VEC_SHUFFLE(r5, r7, F32_SHUFFLE_CONST(1, 0, 1, 0));
      r2 = F32_VEC_SHUFFLE(r5, r7, F32_SHUFFLE_CONST(3, 2, 3, 2));
      r3 = F32_VEC_SHUFFLE(r6, r8, F32_SHUFFLE_CONST(1, 0, 1, 0));
      r4 = F32_VEC_SHUFFLE(r6, r8, F32_SHUFFLE_CONST(3, 2, 3, 2));
		
      *r1_ptr++ = F32_VEC_SHUFFLE(r1, r3, F32_SHUFFLE_CONST(2, 0, 2, 0));
      *r2_ptr++ = F32_VEC_SHUFFLE(r2, r4, F32_SHUFFLE_CONST(2, 0, 2, 0));
      *r3_ptr++ = F32_VEC_SHUFFLE(r1, r3, F32_SHUFFLE_CONST(3, 1, 3, 1));
      *r4_ptr++ = F32_VEC_SHUFFLE(r2, r4, F32_SHUFFLE_CONST(3, 1, 3, 1));
		
      i1 = F32_VEC_SHUFFLE(i5, i7, F32_SHUFFLE_CONST(1, 0, 1, 0));
      i2 = F32_VEC_SHUFFLE(i5, i7, F32_SHUFFLE_CONST(3, 2, 3, 2));
      i3 = F32_VEC_SHUFFLE(i6, i8, F32_SHUFFLE_CONST(1, 0, 1, 0));
      i4 = F32_VEC_SHUFFLE(i6, i8, F32_SHUFFLE_CONST(3, 2, 3, 2));
		
      *i1_ptr++ = F32_VEC_SHUFFLE(i1, i3, F32_SHUFFLE_CONST(2, 0, 2, 0));
      *i2_ptr++ = F32_VEC_SHUFFLE(i2, i4, F32_SHUFFLE_CONST(2, 0, 2, 0));
      *i3_ptr++ = F32_VEC_SHUFFLE(i1, i3, F32_SHUFFLE_CONST(3, 1, 3, 1));
      *i4_ptr++ = F32_VEC_SHUFFLE(i2, i4, F32_SHUFFLE_CONST(3, 1, 3, 1));
    }
}

void pass_3_reorder_simd(FFT_Split *input, t_uint length)
{
  t_uint offset = length >> 5;
  t_uint outer_loop = length >> 6;
  t_uint i, j;
	
  vFloat r1, r2, r3, r4, r5;
  vFloat i1, i2, i3, i4, i5;
  vFloat twiddle_c = {1, (t_sample) SQRT_2_2, 0, (t_sample) -SQRT_2_2};
  vFloat twiddle_s = {0, (t_sample)-SQRT_2_2, -1, (t_sample) -SQRT_2_2};
	
  vFloat *r1_ptr = (vFloat *)input->realp;
  vFloat *i1_ptr = (vFloat *)input->imagp;
  vFloat *r2_ptr = r1_ptr + offset;
  vFloat *i2_ptr = i1_ptr + offset;
	
  for (i = 0, j = 0; i < length >> 1; i += 8)
    {
      // Get input
		
      r1 = *r1_ptr;
      i1 = *i1_ptr;
      r3 = *(r1_ptr + 1);
      i3 = *(i1_ptr + 1);
      r2 = *r2_ptr;
      i2 = *i2_ptr;
      r4 = *(r2_ptr + 1);
      i4 = *(i2_ptr + 1);
		
      // Multiply by twiddle
		
      r5 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r2, twiddle_c),
			  F32_VEC_MUL_OP(i2, twiddle_s));
      i5 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r2, twiddle_s),
			  F32_VEC_MUL_OP(i2, twiddle_c));
		
      // Store output (swapping as necessary)
		
      *r1_ptr = F32_VEC_ADD_OP(r1, r5);
      *i1_ptr = F32_VEC_ADD_OP(i1, i5);
		
      *(r1_ptr + 1) = F32_VEC_SUB_OP(r1, r5);
      *(i1_ptr + 1) = F32_VEC_SUB_OP(i1, i5);
		
      // Multiply by twiddle
		
      r5 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r4, twiddle_c),
			  F32_VEC_MUL_OP(i4, twiddle_s));
      i5 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r4, twiddle_s),
			  F32_VEC_MUL_OP(i4, twiddle_c));
		
      // Store output (swapping as necessary)
		
      *r2_ptr = F32_VEC_ADD_OP(r3, r5);
      *i2_ptr = F32_VEC_ADD_OP(i3, i5);
		
      *(r2_ptr + 1) = F32_VEC_SUB_OP(r3, r5);
      *(i2_ptr + 1) = F32_VEC_SUB_OP(i3, i5);
		
      r1_ptr += 2;
      r2_ptr += 2;
      i1_ptr += 2;
      i2_ptr += 2;
		
      if (!(++j % outer_loop))
	{
	  r1_ptr += offset;
	  r2_ptr += offset;
	  i1_ptr += offset;
	  i2_ptr += offset;
	}
    }
}

void pass_trig_table_reorder_simd(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass)
{
  t_uint size = 2 << pass;
  t_uint incr = size >> 3;
  t_uint loop = size;
  t_uint offset = length >> (pass + 3);
  t_uint outer_loop = ((length >> 1) / size) / ((t_uint) 1 << pass);
  t_uint i, j;
	
  vFloat r1, r2, r3, r4, r5;
  vFloat i1, i2, i3, i4, i5;
  vFloat twiddle_c, twiddle_s;
	
  vFloat *r1_ptr = (vFloat *)input->realp;
  vFloat *i1_ptr = (vFloat *)input->imagp;
  vFloat *r2_ptr = r1_ptr + offset;
  vFloat *i2_ptr = i1_ptr + offset;
	
  for (j = 0, i = 0; i < length >> 1; loop += size)
    {
      vFloat *table_r = (vFloat *)setup->tables[pass - PASS_TRIG_OFFSET].realp;
      vFloat *table_i = (vFloat *)setup->tables[pass - PASS_TRIG_OFFSET].imagp;
		
      for (; i < loop; i += 8)
	{
	  // Get input
			
	  r1 = *r1_ptr;
	  i1 = *i1_ptr;
	  r2 = *r2_ptr;
	  i2 = *i2_ptr;
			
	  // Get Twiddle
			
	  twiddle_c = *table_r++;
	  twiddle_s = *table_i++;
			
	  // Multiply by twiddle
			
	  r5 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r2, twiddle_c),
			      F32_VEC_MUL_OP(i2, twiddle_s));
	  i5 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r2, twiddle_s),
			      F32_VEC_MUL_OP(i2, twiddle_c));
			
	  // Get input
			
	  r3 = *(r1_ptr + incr);
	  i3 = *(i1_ptr + incr);
	  r4 = *(r2_ptr + incr);
	  i4 = *(i2_ptr + incr);
			
	  // Store output (swapping as necessary)
			
	  *r1_ptr = F32_VEC_ADD_OP(r1, r5);
	  *i1_ptr = F32_VEC_ADD_OP(i1, i5);
			
	  *(r1_ptr++ + incr) = F32_VEC_SUB_OP(r1, r5);
	  *(i1_ptr++ + incr) = F32_VEC_SUB_OP(i1, i5);
			
	  // Multiply by twiddle
			
	  r5 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r4, twiddle_c),
			      F32_VEC_MUL_OP(i4, twiddle_s));
	  i5 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r4, twiddle_s),
			      F32_VEC_MUL_OP(i4, twiddle_c));
			
	  // Store output (swapping as necessary)
			
	  *r2_ptr = F32_VEC_ADD_OP(r3, r5);
	  *i2_ptr = F32_VEC_ADD_OP(i3, i5);
			
	  *(r2_ptr++ + incr) = F32_VEC_SUB_OP(r3, r5);
	  *(i2_ptr++ + incr) = F32_VEC_SUB_OP(i3, i5);
	}
		
      r1_ptr += incr;
      r2_ptr += incr;
      i1_ptr += incr;
      i2_ptr += incr;
		
      if (!(++j % outer_loop))
	{
	  r1_ptr += offset;
	  r2_ptr += offset;
	  i1_ptr += offset;
	  i2_ptr += offset;
	}
    }
}

void pass_trig_table_simd(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass)
{

  /* post("vFloat"); */
  
  t_uint size = 2 << pass;
  t_uint incr = size >> 3;
  t_uint loop = size;
  t_uint i;
	
  vFloat r1, r2, r3, i1, i2, i3;
  vFloat twiddle_c, twiddle_s;
	
  vFloat *r1_ptr = (vFloat *)input->realp;
  vFloat *i1_ptr = (vFloat *)input->imagp;
  vFloat *r2_ptr = r1_ptr + incr;
  vFloat *i2_ptr = i1_ptr + incr;
	
  for (i = 0; i < length; loop += size)
    {
      vFloat *table_r = (vFloat *)setup->tables[pass - PASS_TRIG_OFFSET].realp;
      vFloat *table_i = (vFloat *)setup->tables[pass - PASS_TRIG_OFFSET].imagp;
		
      for (; i < loop; i += 8)
	{
	  // Get input
			
	  r1 = *r1_ptr;
	  i1 = *i1_ptr;
	  r2 = *r2_ptr;
	  i2 = *i2_ptr;
			
	  // Get Twiddle
			
	  twiddle_c = *table_r++;
	  twiddle_s = *table_i++;
			
	  // Multiply by twiddle
			
	  r3 = F32_VEC_SUB_OP(F32_VEC_MUL_OP(r2, twiddle_c),
			      F32_VEC_MUL_OP(i2, twiddle_s));
	  i3 = F32_VEC_ADD_OP(F32_VEC_MUL_OP(r2, twiddle_s),
			      F32_VEC_MUL_OP(i2, twiddle_c));
			
	  // Store output (same pos as inputs)
			
	  *r1_ptr++ = F32_VEC_ADD_OP(r1, r3);
	  *i1_ptr++ = F32_VEC_ADD_OP(i1, i3);
			
	  *r2_ptr++ = F32_VEC_SUB_OP(r1, r3);
	  *i2_ptr++ = F32_VEC_SUB_OP(i1, i3);
	}
		
      r1_ptr += incr;
      r2_ptr += incr;
      i1_ptr += incr;
      i2_ptr += incr;
    }
}
#else

void pass_1_2_reorder_simd(FFT_Split *input, t_uint length)
{
  t_uint i;

  vDouble *r1_ptr = (vDouble *)input->realp;
  vDouble *r2_ptr = r1_ptr + (length >> 3);
  vDouble *r3_ptr = r2_ptr + (length >> 3);
  vDouble *r4_ptr = r3_ptr + (length >> 3);

  vDouble *i1_ptr = (vDouble *)input->imagp;
  vDouble *i2_ptr = i1_ptr + (length >> 3);
  vDouble *i3_ptr = i2_ptr + (length >> 3);
  vDouble *i4_ptr = i3_ptr + (length >> 3);

  vDouble r1, r2, r3, r4, r5, r6, r7, r8;
  vDouble i1, i2, i3, i4, i5, i6, i7, i8;
  vDouble t1, t2, t3, t4;

  for(i=0; i<length>>4; i++) {
    t1 = *r1_ptr;
    t2 = *(r1_ptr + 1);
    t3 = *(r2_ptr);
    t4 = *(r2_ptr + 1);
    r5 = *(r3_ptr);
    r6 = *(r3_ptr + 1);
    r7 = *(r4_ptr);
    r8 = *(r4_ptr + 1);

    r1 = F64_VEC_ADD_OP(t1, r5);
    r2 = F64_VEC_ADD_OP(t2, r6);
    r3 = F64_VEC_ADD_OP(t3, r7);
    r4 = F64_VEC_ADD_OP(t4, r8);
    r5 = F64_VEC_SUB_OP(t1, r5);
    r6 = F64_VEC_SUB_OP(t2, r6);
    r7 = F64_VEC_SUB_OP(t3, r7);
    r8 = F64_VEC_SUB_OP(t4, r8);

    t1 = *i1_ptr;
    t2 = *(i1_ptr + 1);
    t3 = *i2_ptr;
    t4 = *(i2_ptr + 1);
    i5 = *i3_ptr;
    i6 = *(i3_ptr + 1);
    i7 = *i4_ptr;
    i8 = *(i4_ptr + 1);

    i1 = F64_VEC_ADD_OP(t1, i5);
    i2 = F64_VEC_ADD_OP(t2, i6);
    i3 = F64_VEC_ADD_OP(t3, i7);
    i4 = F64_VEC_ADD_OP(t4, i8);
    i5 = F64_VEC_SUB_OP(t1, i5);
    i6 = F64_VEC_SUB_OP(t2, i6);
    i7 = F64_VEC_SUB_OP(t3, i7);
    i8 = F64_VEC_SUB_OP(t4, i8);

    t1 = F64_VEC_ADD_OP(r1, r3);
    t2 = F64_VEC_ADD_OP(r2, r4);
    t3 = F64_VEC_SUB_OP(r1, r3);
    t4 = F64_VEC_SUB_OP(r2, r4);
    r1 = F64_VEC_ADD_OP(r5, i7);
    r2 = F64_VEC_ADD_OP(r6, i8);
    r3 = F64_VEC_SUB_OP(r5, i7);
    r4 = F64_VEC_SUB_OP(r6, i8);

    *r1_ptr++ = F64_VEC_SHUFFLE(t1, r1, F64_SHUFFLE_CONST(0, 0));
    *r1_ptr++ = F64_VEC_SHUFFLE(t3, r3, F64_SHUFFLE_CONST(0, 0));
    *r2_ptr++ = F64_VEC_SHUFFLE(t2, r2, F64_SHUFFLE_CONST(0, 0));
    *r2_ptr++ = F64_VEC_SHUFFLE(t4, r4, F64_SHUFFLE_CONST(0, 0));
    *r3_ptr++ = F64_VEC_SHUFFLE(t1, r1, F64_SHUFFLE_CONST(1, 1));
    *r3_ptr++ = F64_VEC_SHUFFLE(t3, r3, F64_SHUFFLE_CONST(1, 1));
    *r4_ptr++ = F64_VEC_SHUFFLE(t2, r2, F64_SHUFFLE_CONST(1, 1));
    *r4_ptr++ = F64_VEC_SHUFFLE(t4, r4, F64_SHUFFLE_CONST(1, 1));

    t1 = F64_VEC_ADD_OP(i1, i3);
    t2 = F64_VEC_ADD_OP(i2, i4);
    t3 = F64_VEC_SUB_OP(i1, i3);
    t4 = F64_VEC_SUB_OP(i2, i4);
    i1 = F64_VEC_SUB_OP(i5, r7);
    i2 = F64_VEC_SUB_OP(i6, r8);
    i3 = F64_VEC_ADD_OP(i5, r7);
    i4 = F64_VEC_ADD_OP(i6, r8);

    *i1_ptr++ = F64_VEC_SHUFFLE(t1, i1, F64_SHUFFLE_CONST(0, 0));
    *i1_ptr++ = F64_VEC_SHUFFLE(t3, i3, F64_SHUFFLE_CONST(0, 0));
    *i2_ptr++ = F64_VEC_SHUFFLE(t2, i2, F64_SHUFFLE_CONST(0, 0));
    *i2_ptr++ = F64_VEC_SHUFFLE(t4, i4, F64_SHUFFLE_CONST(0, 0));
    *i3_ptr++ = F64_VEC_SHUFFLE(t1, i1, F64_SHUFFLE_CONST(1, 1));
    *i3_ptr++ = F64_VEC_SHUFFLE(t3, i3, F64_SHUFFLE_CONST(1, 1));
    *i4_ptr++ = F64_VEC_SHUFFLE(t2, i2, F64_SHUFFLE_CONST(1, 1));
    *i4_ptr++ = F64_VEC_SHUFFLE(t4, i4, F64_SHUFFLE_CONST(1, 1));
  }
}

void pass_3_reorder_simd(FFT_Split *input, t_uint length)
{
  t_uint offset = length >> 4;
  t_uint outer_loop = length >> 6;
  t_uint i, j;

  vDouble r1, r2, r3, r4, r5, r6, r7, r8, r9, r10;
  vDouble i1, i2, i3, i4, i5, i6, i7, i8, i9, i10;

  vDouble twiddle_c1 = {1, SQRT_2_2};
  vDouble twiddle_s1 = {0, -SQRT_2_2};
  vDouble twiddle_c2 = {0, -SQRT_2_2};
  vDouble twiddle_s2 = {-1, -SQRT_2_2};

  vDouble *r1_ptr = (vDouble *)input->realp;
  vDouble *i1_ptr = (vDouble *)input->imagp;
  vDouble *r2_ptr = r1_ptr + offset;
  vDouble *i2_ptr = i1_ptr + offset;

  for(i = 0, j = 0; i<length >> 1; i+= 8) {
    r1 = *(r1_ptr);
    r5 = *(r1_ptr + 1);
    r2 = *r2_ptr;
    r6 = *(r2_ptr + 1);
    i1 = *i1_ptr;
    i5 = *(i1_ptr + 1);
    i2 = *i2_ptr;
    i6 = *(i2_ptr + 1);

    r9 = F64_VEC_SUB_OP(F64_VEC_MUL_OP(r2, twiddle_c1),
			F64_VEC_MUL_OP(i2, twiddle_s1));
    i9 = F64_VEC_ADD_OP(F64_VEC_MUL_OP(r2, twiddle_s1),
			F64_VEC_MUL_OP(i2, twiddle_c1));
    r10 = F64_VEC_SUB_OP(F64_VEC_MUL_OP(r6, twiddle_c2),
			 F64_VEC_MUL_OP(i6, twiddle_s2));
    i10 = F64_VEC_ADD_OP(F64_VEC_MUL_OP(r6, twiddle_s2),
			 F64_VEC_MUL_OP(i6, twiddle_c2));

    r3 = *(r1_ptr + 2);
    r7 = *(r1_ptr + 3);
    r4 = *(r2_ptr + 2);
    r8 = *(r2_ptr + 3);
    i3 = *(i1_ptr + 2);
    i7 = *(i1_ptr + 3);
    i4 = *(i2_ptr + 2);
    i8 = *(i2_ptr + 3);

    *r1_ptr = F64_VEC_ADD_OP(r1, r9);
    *(r1_ptr + 1) = F64_VEC_ADD_OP(r5, r10);
    *i1_ptr = F64_VEC_ADD_OP(i1, i9);
    *(i1_ptr + 1) = F64_VEC_ADD_OP(i5, i10);

    *(r1_ptr++ + 2) = F64_VEC_SUB_OP(r1, r9);
    *(r1_ptr++ + 2) = F64_VEC_SUB_OP(r5, r10);
    *(i1_ptr++ + 2) = F64_VEC_SUB_OP(i1, i9);
    *(i1_ptr++ + 2) = F64_VEC_SUB_OP(i5, i10);

    r9 = F64_VEC_SUB_OP(F64_VEC_MUL_OP(r4, twiddle_c1),
			F64_VEC_MUL_OP(i4, twiddle_s1));
    i9 = F64_VEC_ADD_OP(F64_VEC_MUL_OP(r4, twiddle_s1),
			F64_VEC_MUL_OP(i4, twiddle_c1));
    r10 = F64_VEC_SUB_OP(F64_VEC_MUL_OP(r8, twiddle_c2),
			 F64_VEC_MUL_OP(i8, twiddle_s2));
    i10 = F64_VEC_ADD_OP(F64_VEC_MUL_OP(r8, twiddle_s2),
			 F64_VEC_MUL_OP(i8, twiddle_c2));

    *r2_ptr = F64_VEC_ADD_OP(r3, r9);
    *(r2_ptr + 1) = F64_VEC_ADD_OP(r7, r10);
    *i2_ptr = F64_VEC_ADD_OP(i3, i9);
    *(i2_ptr + 1) = F64_VEC_ADD_OP(i7, i10);

    *(r2_ptr++ + 2) = F64_VEC_SUB_OP(r3, r9);
    *(r2_ptr++ + 2) = F64_VEC_SUB_OP(r7, r10);
    *(i2_ptr++ + 2) = F64_VEC_SUB_OP(i3, i9);
    *(i2_ptr++ + 2) = F64_VEC_SUB_OP(i7, i10);

    r1_ptr += 2;
    r2_ptr += 2;
    i1_ptr += 2;
    i2_ptr += 2;

    if(!(++j % outer_loop)) {
      r1_ptr += offset;
      r2_ptr += offset;
      i1_ptr += offset;
      i2_ptr += offset;
    }
  }
}

void pass_trig_table_reorder_simd(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass)
{
  t_uint size = 2 << pass;
  t_uint incr = size >> 2;
  t_uint loop = size;
  t_uint offset = length >> (pass + 2);
  t_uint outer_loop = ((length >> 1) / size) / ((t_uint) 1 << pass);
  t_uint i, j;

  vDouble r1, r2, r3, r4, r5;
  vDouble i1, i2, i3, i4, i5;
  vDouble twiddle_c, twiddle_s;

  vDouble *r1_ptr = (vDouble *)input->realp;
  vDouble *i1_ptr = (vDouble *)input->imagp;
  vDouble *r2_ptr = r1_ptr + offset;
  vDouble *i2_ptr = i1_ptr + offset;

  for(j = 0, i = 0; i < (length >> 1); loop += size) {
    vDouble *t_ptr_r = (vDouble *)setup->tables[pass - PASS_TRIG_OFFSET].realp;
    vDouble *t_ptr_i = (vDouble *)setup->tables[pass - PASS_TRIG_OFFSET].imagp;

    for(; i<loop; i+= 4) {
      r1 = *r1_ptr;
      i1 = *i1_ptr;
      r2 = *r2_ptr;
      i2 = *i2_ptr;

      twiddle_c = *t_ptr_r++;
      twiddle_s = *t_ptr_i++;

      r5 = F64_VEC_SUB_OP(F64_VEC_MUL_OP(r2, twiddle_c),
			  F64_VEC_MUL_OP(i2, twiddle_s));
      i5 = F64_VEC_ADD_OP(F64_VEC_MUL_OP(r2, twiddle_s),
			  F64_VEC_MUL_OP(i2, twiddle_c));

      r3 = *(r1_ptr + incr);
      i3 = *(i1_ptr + incr);
      r4 = *(r2_ptr + incr);
      i4 = *(i2_ptr + incr);

      *r1_ptr = F64_VEC_ADD_OP(r1, r5);
      *i1_ptr = F64_VEC_ADD_OP(i1, i5);

      *(r1_ptr++ + incr) = F64_VEC_SUB_OP(r1, r5);
      *(i1_ptr++ + incr) = F64_VEC_SUB_OP(i1, i5);

      r5 = F64_VEC_SUB_OP(F64_VEC_MUL_OP(r4, twiddle_c),
			  F64_VEC_MUL_OP(i4, twiddle_s));
      i5 = F64_VEC_ADD_OP(F64_VEC_MUL_OP(r4, twiddle_s),
			  F64_VEC_MUL_OP(i4, twiddle_c));

      *r2_ptr = F64_VEC_ADD_OP(r3, r5);
      *i2_ptr = F64_VEC_ADD_OP(i3, i5);

      *(r2_ptr++ + incr) = F64_VEC_SUB_OP(r3, r5);
      *(i2_ptr++ + incr) = F64_VEC_SUB_OP(i3, i5);
    }

    r1_ptr += incr;
    r2_ptr += incr;
    i1_ptr += incr;
    i2_ptr += incr;

    if(!(++j % outer_loop)) {
      r1_ptr += offset;
      r2_ptr += offset;
      i1_ptr += offset;
      i2_ptr += offset;
    }
  }
}

void pass_trig_table_simd(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass)
{

  /* post("vDouble"); */
  
  t_uint size = 2 << pass;
  t_uint incr = size >> 2;
  t_uint loop = size;

  vDouble r1, r2, i1, i2, r3, i3, twiddle_c, twiddle_s;

  vDouble *r1_ptr = (vDouble *)input->realp;
  vDouble *i1_ptr = (vDouble *)input->imagp;
  vDouble *r2_ptr = r1_ptr + incr;
  vDouble *i2_ptr = i1_ptr + incr;

  t_uint i;

  for(i=0; i<length; loop += size) {
    vDouble *t_ptr_r = (vDouble *)setup->tables[pass - PASS_TRIG_OFFSET].realp;
    vDouble *t_ptr_i = (vDouble *)setup->tables[pass - PASS_TRIG_OFFSET].imagp;

    for(; i < loop; i += 4) {
      r1 = *r1_ptr;
      i1 = *i1_ptr;
      r2 = *r2_ptr;
      i2 = *i2_ptr;

      twiddle_c = *t_ptr_r++;
      twiddle_s = *t_ptr_i++;

      r3 = F64_VEC_SUB_OP(F64_VEC_MUL_OP(r2, twiddle_c),
			  F64_VEC_MUL_OP(i2, twiddle_s));
      i3 = F64_VEC_ADD_OP(F64_VEC_MUL_OP(r2, twiddle_s),
			  F64_VEC_MUL_OP(i2, twiddle_c));

      *r1_ptr++ = F64_VEC_ADD_OP(r1, r3);
      *i1_ptr++ = F64_VEC_ADD_OP(i1, i3);

      *r2_ptr++ = F64_VEC_SUB_OP(r1, r3);
      *i2_ptr++ = F64_VEC_SUB_OP(i1, i3);
    }

    r1_ptr += incr;
    r2_ptr += incr;
    i1_ptr += incr;
    i2_ptr += incr;
  }
}
#endif /* PD_SAMPLE_PRECISION */


void pass_real_trig_table32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2, long ifft)
{
  t_uint length = (t_uint) 1 << (fft_log2 - 1);
  t_uint length_m1 = length - 1;

  t_float32 r1, r2, r3, r4, i1, i2, i3, i4, t1, t2;
  t_float32 twiddle_c1, twiddle_s1;

  t_float32 *r1_ptr = input->realp;
  t_float32 *i1_ptr = input->imagp;

  t_float32 *r2_ptr = r1_ptr + length_m1;
  t_float32 *i2_ptr = i1_ptr + length_m1;
  t_float32 *tr1_ptr = setup->tables[fft_log2 - FFTLOG2_TRIG_OFFSET].realp;
  t_float32 *ti1_ptr = setup->tables[fft_log2 - FFTLOG2_TRIG_OFFSET].imagp;

  t_float32 flip = 1.;

  t_uint i;

  if(ifft)
    flip = -1.;

  tr1_ptr++;
  ti1_ptr++;

  r1 = *r1_ptr;
  i1 = *i1_ptr;

  t1 = r1 + i1;
  t2 = r1 - i1;

  if(!ifft) {
    t1 *= 2.;
    t2 *= 2.;
  }

  *r1_ptr++ = t1;
  *i1_ptr++ = t2;

  for(i=0; i < (length >> 1); i++) {
    twiddle_c1 = flip * *tr1_ptr++;
    twiddle_s1 = *ti1_ptr++;

    r1 = *r1_ptr;
    i1 = *i1_ptr;
    r2 = *r2_ptr;
    i2 = *i2_ptr;

    r3 = r1 + r2;
    i3 = i1 + i2;
    r4 = r1 - r2;
    i4 = i1 - i2;

    t1 = (twiddle_c1 * i3) + (twiddle_s1 * r4);
    t2 = (twiddle_c1 * -r4) + (twiddle_s1 * i3);

    *r1_ptr++ = r3 + t1;
    *i1_ptr++ = t2 + i4;

    *r2_ptr-- = r3 - t1;
    *i2_ptr-- = t2 - i4;
  }
}

void do_small_real_fft32(FFT_Split32 *input, t_uint fft_log2, long ifft)
{
  t_float32 r1, r2, r3, r4, i1, i2;
  t_float32 scale = 2.;

  t_float32 *realp = input->realp;
  t_float32 *imagp = input->imagp;

  if(fft_log2 < 1)
    return;

  if(ifft)
    scale = 1.;

  switch(fft_log2)
    {
    case 1:
      r1 = realp[0];
      r2 = imagp[0];

      realp[0] = (r1 + r2) * scale;
      imagp[0] = (r1 - r2) * scale;

      break;

    case 2:
      if(!ifft) {
	r3 = realp[0];
	r4 = realp[1];
	i1 = imagp[0];
	i2 = imagp[1];

	r1 = r3 + r4;
	r2 = r3 - r4;
	r3 = i1 + i2;
	r4 = i1 - i2;

	realp[0] = (r1 + r3) * 2;
	realp[1] = r2 * 2;
	imagp[0] = (r1 - r3) * 2;
	imagp[1] = -r4 * 2;
      }
      else {
	i1 = realp[0];
	r2 = realp[1] + realp[1];
	i2 = imagp[0];
	r4 = -imagp[1] - imagp[1];

	r1 = i1 + i2;
	r3 = i1 - i2;

	realp[0] = r1 + r2;
	realp[1] = r1 - r2;
	imagp[0] = r3 + r4;
	imagp[1] = r3 - r4;
      }
      break;
    }
}

void pass_trig_table32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint length, t_uint pass)
{
  t_uint size = 2 << pass;
  t_uint incr = size >> 1;
  t_uint loop = size;

  t_uint i;

  double r1, r2, i1, i2, r_0, i_0;
  double twiddle_c, twiddle_s;

  t_float32 *r1_ptr = (t_float32 *)input->realp;
  t_float32 *i1_ptr = (t_float32 *)input->imagp;
  t_float32 *r2_ptr = (t_float32 *)r1_ptr + incr;
  t_float32 *i2_ptr = (t_float32 *)i1_ptr + incr;

  for(i = 0; i < length; loop += size) {
    t_float32 *tr_ptr = setup->tables[pass - PASS_TRIG_OFFSET].realp;
    t_float32 *ti_ptr = setup->tables[pass - PASS_TRIG_OFFSET].imagp;

    for(; i < loop; i += 2) {
      twiddle_c = (double)*tr_ptr++;
      twiddle_s = (double)*ti_ptr++;

      r1 = *r1_ptr;
      i1 = *i1_ptr;
      r2 = *r2_ptr;
      i2 = *i2_ptr;

      r_0 = (r2 * twiddle_c) - (i2 * twiddle_s);
      i_0 = (r2 * twiddle_s) + (i2 * twiddle_c);

      *r1_ptr++ = (t_sample) (r1 + r_0);
      *i1_ptr++ = (t_sample) (i1 + i_0);

      *r2_ptr++ = (t_sample) (r1 - r_0);
      *i2_ptr++ = (t_sample) (i1 - i_0);
    }
    
    r1_ptr += incr;
    r2_ptr += incr;
    i1_ptr += incr;
    i2_ptr += incr;
  }
}

void pass_trig_table_reorder32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint length, t_uint pass)
{
  t_uint size = 2 << pass;
  t_uint incr = size >> 1;
  t_uint loop = size;
  t_uint offset = length >> (pass + 1);
  t_uint outer_loop = ((length >> 1) / size ) / ((t_uint) 1 << pass);
  t_uint i, j;
    
  t_float32 r1, r2, r3, r4, r5;
  t_float32 i1, i2, i3, i4, i5;
  t_float32 twiddle_c, twiddle_s;

  t_float32 *r1_ptr = input->realp;
  t_float32 *i1_ptr = input->imagp;
  t_float32 *r2_ptr = r1_ptr + incr;
  t_float32 *i2_ptr = i1_ptr + incr;
  t_float32 *r3_ptr = r1_ptr + offset;
  t_float32 *i3_ptr = i1_ptr + offset;
  t_float32 *r4_ptr = r3_ptr + incr;
  t_float32 *i4_ptr = i3_ptr + incr;

  for(j = 0, i = 0; i < (length >> 1); loop += size) {
    t_float32 *t_ptr_r = setup->tables[pass - PASS_TRIG_OFFSET].realp;
    t_float32 *t_ptr_i = setup->tables[pass - PASS_TRIG_OFFSET].imagp;

    for(; i < loop; i += 2) {
      r1 = *r1_ptr;
      i1 = *i1_ptr;
      r2 = *r3_ptr;
      i2 = *i3_ptr;

      twiddle_c = *t_ptr_r++;
      twiddle_s = *t_ptr_i++;

      r5 = (r2 * twiddle_c) - (i2 * twiddle_s);
      i5 = (r2 * twiddle_s) + (i2 * twiddle_c);

      r3 = *r2_ptr;
      i3 = *i2_ptr;
      r4 = *r4_ptr;
      i4 = *i4_ptr;

      *r1_ptr++ = r1 + r5;
      *i1_ptr++ = i1 + i5;

      *r2_ptr++ = r1 - r5;
      *i2_ptr++ = i1 - i5;

      r5 = (r4 * twiddle_c) - (i4 * twiddle_s);
      i5 = (r4 * twiddle_s) + (i4 * twiddle_c);

      *r3_ptr++ = r3 + r5;
      *i3_ptr++ = i3 + i5;

      *r4_ptr++ = r3 - r5;
      *i4_ptr++ = i3 - i5;
    }

    r1_ptr += incr;
    r2_ptr += incr;
    r3_ptr += incr;
    r4_ptr += incr;
    i1_ptr += incr;
    i2_ptr += incr;
    i3_ptr += incr;
    i4_ptr += incr;

    if(!(++j % outer_loop)) {
      r1_ptr += offset;
      r2_ptr += offset;
      r3_ptr += offset;
      r4_ptr += offset;
      i1_ptr += offset;
      i2_ptr += offset;
      i3_ptr += offset;
      i4_ptr += offset;
    }
  }
}

void pass_3_reorder32(FFT_Split32 *input, t_uint length, t_uint fft_log2)
{
  t_uint offset = (t_uint) 1 << (fft_log2 - 3);
  t_uint i, j, loop;

  t_float32 sqrt_2_2 = (t_float32) SQRT_2_2;
  t_float32 r1, r2, r3, r4, r5, r6, r7, r8;
  t_float32 r9, r10, r11, r12, r13, r14, r15, r16;
  t_float32 i1, i2, i3, i4, i5, i6, i7, i8;
  t_float32 i9, i10, i11, i12, i13, i14, i15, i16;
  t_float32 t1, t2, t3, t4;

  t_float32 *r1_ptr = input->realp;
  t_float32 *i1_ptr = input->imagp;

  t_float32 *r2_ptr = r1_ptr + offset;
  t_float32 *i2_ptr = i1_ptr + offset;

  for(j = 0, loop = length >> 6; j < 4; j++) {
    for(i=0; i<loop; i++) {
      r1  = *(r1_ptr + 0);
      r2  = *(r1_ptr + 1);
      r3  = *(r1_ptr + 2);
      r4  = *(r1_ptr + 3);

      r9  = *(r1_ptr + 4);
      r10 = *(r1_ptr + 5);
      r11 = *(r1_ptr + 6);
      r12 = *(r1_ptr + 7);

      r5  = *(r2_ptr + 0);
      r6  = *(r2_ptr + 1);
      r7  = *(r2_ptr + 2);
      r8  = *(r2_ptr + 3);

      r13 = *(r2_ptr + 4);
      r14 = *(r2_ptr + 5);
      r15 = *(r2_ptr + 6);
      r16 = *(r2_ptr + 7);

      i1  = *(i1_ptr + 0);
      i2  = *(i1_ptr + 1);
      i3  = *(i1_ptr + 2);
      i4  = *(i1_ptr + 3);

      i9  = *(i1_ptr + 4);
      i10 = *(i1_ptr + 5);
      i11 = *(i1_ptr + 6);
      i12 = *(i1_ptr + 7);

      i5  = *(i2_ptr + 0);
      i6  = *(i2_ptr + 1);
      i7  = *(i2_ptr + 2);
      i8  = *(i2_ptr + 3);

      i13 = *(i2_ptr + 4);
      i14 = *(i2_ptr + 5);
      i15 = *(i2_ptr + 6);
      i16 = *(i2_ptr + 7);

      t1 = sqrt_2_2 * (r6 + i6);
      t2 = sqrt_2_2 * (i8 - r8);
      t3 = sqrt_2_2 * (r14 + i14);
      t4 = sqrt_2_2 * (i16 - r16);

      *r1_ptr++ = r1 + r5;
      *r1_ptr++ = r2 + t1;
      *r1_ptr++ = r3 + i7;
      *r1_ptr++ = r4 + t2;
      *r1_ptr++ = r1 - r5;
      *r1_ptr++ = r2 - t1;
      *r1_ptr++ = r3 - i7;
      *r1_ptr++ = r4 - t2;

      *r2_ptr++ = r9  + r13;
      *r2_ptr++ = r10 + t3;
      *r2_ptr++ = r11 + i15;
      *r2_ptr++ = r12 + t4;
      *r2_ptr++ = r9  - r13;
      *r2_ptr++ = r10 - t3;
      *r2_ptr++ = r11 - i15;
      *r2_ptr++ = r12 - t4;

      t1 = sqrt_2_2 * (i6 - r6);
      t2 = -sqrt_2_2 * (r8 + i8);
      t3 = sqrt_2_2 * (i14 - r14);
      t4 = -sqrt_2_2 * (r16 + i16);

      *i1_ptr++ = i1 + i5;
      *i1_ptr++ = i2 + t1;
      *i1_ptr++ = i3 - r7;
      *i1_ptr++ = i4 + t2;
      *i1_ptr++ = i1 - i5;
      *i1_ptr++ = i2 - t1;
      *i1_ptr++ = i3 + r7;
      *i1_ptr++ = i4 - t2;

      *i2_ptr++ = i9  + i13;
      *i2_ptr++ = i10 + t3;
      *i2_ptr++ = i11 - r15;
      *i2_ptr++ = i12 + t4;
      *i2_ptr++ = i9  - i13;
      *i2_ptr++ = i10 - t3;
      *i2_ptr++ = i11 + r15;
      *i2_ptr++ = i12 - t4;
    }
    r1_ptr += offset;
    i1_ptr += offset;
    r2_ptr += offset;
    i2_ptr += offset;
  }
}

void pass_1_2_reorder32(FFT_Split32 *input, t_uint length, t_uint fft_log2)
{
  t_uint i, j, loop;
  t_uint offset = length >> 1;

  t_float32 *r1_ptr = input->realp;
  t_float32 *i1_ptr = input->imagp;
  t_float32 *r2_ptr = r1_ptr + offset;
  t_float32 *i2_ptr = i1_ptr + offset;

  t_float32 r1, r2, r3, r4, r5, r6, r7, r8;
  t_float32 i1, i2, i3, i4, i5, i6, i7, i8;

  for(i=0; i< (length >> 2); i++) {
    r1 = *r1_ptr;
    r2 = *r2_ptr;
    r3 = *(r1_ptr + 1);
    r4 = *(r2_ptr + 1);

    *r1_ptr++ = r1 + r2;
    *r1_ptr++ = r1 - r2;
    *r2_ptr++ = r3 + r4;
    *r2_ptr++ = r3 - r4;

    i1 = *i1_ptr;
    i2 = *i2_ptr;
    i3 = *(i1_ptr + 1);
    i4 = *(i2_ptr + 1);

    *i1_ptr++ = i1 + i2;
    *i1_ptr++ = i1 - i2;
    *i2_ptr++ = i3 + i4;
    *i2_ptr++ = i3 - i4;
  }

  offset >>= 1;

  r1_ptr = input->realp;
  i1_ptr = input->imagp;
  r2_ptr = r1_ptr + offset;
  i2_ptr = i1_ptr + offset;

  for(j=0, loop = (length >> 4); j < 2; j++) {
    for(i=0;i<loop;i++) {
      r1 = *(r1_ptr + 0);
      r2 = *(r1_ptr + 1);
      r3 = *(r2_ptr + 0);
      r4 = *(r2_ptr + 1);
      r5 = *(r1_ptr + 2);
      r6 = *(r1_ptr + 3);
      r7 = *(r2_ptr + 2);
      r8 = *(r2_ptr + 3);

      i1 = *(i1_ptr + 0);
      i2 = *(i1_ptr + 1);
      i3 = *(i2_ptr + 0);
      i4 = *(i2_ptr + 1);
      i5 = *(i1_ptr + 2);
      i6 = *(i1_ptr + 3);
      i7 = *(i2_ptr + 2);
      i8 = *(i2_ptr + 3);

      *r1_ptr++ = r1 + r3;
      *r1_ptr++ = r2 + i4;
      *r1_ptr++ = r1 - r3;
      *r1_ptr++ = r2 - i4;
      *r2_ptr++ = r5 + r7;
      *r2_ptr++ = r6 + i8;
      *r2_ptr++ = r5 - r7;
      *r2_ptr++ = r6 - i8;

      *i1_ptr++ = i1 + i3;
      *i1_ptr++ = i2 - r4;
      *i1_ptr++ = i1 - i3;
      *i1_ptr++ = i2 + r4;
      *i2_ptr++ = i5 + i7;
      *i2_ptr++ = i6 - r8;
      *i2_ptr++ = i5 - i7;
      *i2_ptr++ = i6 + r8;      
    }

    r1_ptr += offset;
    i1_ptr += offset;
    r2_ptr += offset;
    i2_ptr += offset;
  }
}

void pass_3_32(FFT_Split32 *input, t_uint length, t_uint fft_log2)
{
  t_uint i;

  t_float32 sqrt_2_2 = (t_float32) SQRT_2_2;
  t_float32 r1, r2, r3, r4, r5, r6, r7, r8;
  t_float32 i1, i2, i3, i4, i5, i6, i7, i8;
  t_float32 t1, t2;

  t_float32 *r1_ptr = input->realp;
  t_float32 *i1_ptr = input->imagp;

  for(i=0; i < (length >> 3); i++) {
    r1 = *(r1_ptr + 0);
    r2 = *(r1_ptr + 1);
    r3 = *(r1_ptr + 2);
    r4 = *(r1_ptr + 3);
    r5 = *(r1_ptr + 4);
    r6 = *(r1_ptr + 5);
    r7 = *(r1_ptr + 6);
    r8 = *(r1_ptr + 7);

    i1 = *(i1_ptr + 0);
    i2 = *(i1_ptr + 1);
    i3 = *(i1_ptr + 2);
    i4 = *(i1_ptr + 3);
    i5 = *(i1_ptr + 4);
    i6 = *(i1_ptr + 5);
    i7 = *(i1_ptr + 6);
    i8 = *(i1_ptr + 7);

    t1 = sqrt_2_2 * (r6 + i6);
    t2 = sqrt_2_2 * (i8 - r8);

    *r1_ptr++ = r1 + r5;
    *r1_ptr++ = r2 + t1;
    *r1_ptr++ = r3 + i7;
    *r1_ptr++ = r4 + t2;
    *r1_ptr++ = r1 - r5;
    *r1_ptr++ = r2 - t1;
    *r1_ptr++ = r3 - i7;
    *r1_ptr++ = r4 - t2;

    t1 = sqrt_2_2 * (i6 - r6);
    t2 = -sqrt_2_2 * (r8 + i8);

    *i1_ptr++ = i1 + i5;
    *i1_ptr++ = i2 + t1;
    *i1_ptr++ = i3 - r7;
    *i1_ptr++ = i4 + t2;
    *i1_ptr++ = i1 - i5;
    *i1_ptr++ = i2 - t1;
    *i1_ptr++ = i3 + r7;
    *i1_ptr++ = i4 - t2;
  }
}

void do_small_fft32(FFT_Split32 *input, t_uint fft_log2)
{
  t_float32 r1, r2, r3, r4, r5, r6, r7, r8;
  t_float32 i1, i2, i3, i4, i5, i6, i7, i8;
  t_float32 t1, t2, t3, t4;

  t_float32 *realp = input->realp;
  t_float32 *imagp = input->imagp;

  if(fft_log2 < 1)
    return;

  switch(fft_log2)
    {
    case 1:
      
      r1 = realp[0];
      r2 = realp[1];
      i1 = imagp[0];
      i2 = imagp[1];

      realp[0] = r1 + r2;
      realp[1] = r1 - r2;
      imagp[0] = i1 + i2;
      imagp[1] = i1 - i2;

      break;

    case 2:

      r5 = realp[0];
      r6 = realp[1];
      r2 = realp[2];
      r4 = realp[3];
      i5 = imagp[0];
      i6 = imagp[1];
      i2 = imagp[2];
      i4 = imagp[3];

      r1 = r5 + r2;
      r2 = r5 - r2;
      r3 = r6 + r4;
      r4 = r6 - r4;
      i1 = i5 + i2;
      i2 = i5 - i2;
      i3 = i6 + i4;
      i4 = i6 - i4;

      realp[0] = r1 + r3;
      realp[1] = r2 + i4;
      realp[2] = r1 - r3;
      realp[3] = r2 - i4;
      imagp[0] = i1 + i3;
      imagp[1] = i2 - r4;
      imagp[2] = i1 - i3;
      imagp[3] = i2 + i4;

      break;

    case 3:

      t1 = realp[0];
      t3 = realp[1];
      t2 = realp[2];
      t4 = realp[3];
      r2 = realp[4];
      r6 = realp[5];
      r4 = realp[6];
      r8 = realp[7];

      r1 = t1 + r2;
      r2 = t1 - r2;
      r3 = t2 + r4;
      r4 = t2 - r4;
      r5 = t3 + r6;
      r6 = t3 - r6;
      r7 = t4 + r8;
      r8 = t4 - r8;

      t1 = imagp[0];
      t3 = imagp[1];
      t2 = imagp[2];
      t4 = imagp[3];
      i2 = imagp[4];
      i6 = imagp[5];
      i4 = imagp[6];
      i8 = imagp[7];

      i1 = t1 + i2;
      i2 = t1 - i2;
      i3 = t2 + i4;
      i4 = t2 - i4;
      i5 = t3 + i6;
      i6 = t3 - i6;
      i7 = t4 + i8;
      i8 = t4 - i8;

      realp[0] = r1 + r3;
      realp[1] = r2 + i4;
      realp[2] = r1 - r3;
      realp[3] = r2 - i4;
      realp[4] = r5 + r7;
      realp[5] = r6 + i8;
      realp[6] = r5 - r7;
      realp[7] = r6 - i8;

      imagp[0] = i1 + i3;
      imagp[1] = i2 - r4;
      imagp[2] = i1 - i3;
      imagp[3] = i2 + r4;
      imagp[4] = i5 + i7;
      imagp[5] = i6 - r8;
      imagp[6] = i5 - i7;
      imagp[7] = i6 + r8;

      pass_3_32(input, (t_uint)8, (t_uint)3);
 
      break;
    }
}

void pass_real_trig_table(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2, long ifft)
{
  t_uint length = (t_uint) 1 << (fft_log2 - 1);
  t_uint length_m1 = length - 1;

  t_sample r1, r2, r3, r4, i1, i2, i3, i4, t1, t2;
  t_sample twiddle_c1, twiddle_s1;

  t_sample *r1_ptr = input->realp;
  t_sample *i1_ptr = input->imagp;

  t_sample *r2_ptr = r1_ptr + length_m1;
  t_sample *i2_ptr = i1_ptr + length_m1;
  t_sample *tr1_ptr = setup->tables[fft_log2 - FFTLOG2_TRIG_OFFSET].realp;
  t_sample *ti1_ptr = setup->tables[fft_log2 - FFTLOG2_TRIG_OFFSET].imagp;

  t_sample flip = 1.;

  t_uint i;

  if(ifft)
    flip = -1.;

  tr1_ptr++;
  ti1_ptr++;

  r1 = *r1_ptr;
  i1 = *i1_ptr;

  t1 = r1 + i1;
  t2 = r1 - i1;

  if(!ifft) {
    t1 *= 2.;
    t2 *= 2.;
  }

  *r1_ptr++ = t1;
  *i1_ptr++ = t2;

  for(i=0; i < (length >> 1); i++) {
    twiddle_c1 = flip * *tr1_ptr++;
    twiddle_s1 = *ti1_ptr++;

    r1 = *r1_ptr;
    i1 = *i1_ptr;
    r2 = *r2_ptr;
    i2 = *i2_ptr;

    r3 = r1 + r2;
    i3 = i1 + i2;
    r4 = r1 - r2;
    i4 = i1 - i2;

    t1 = (twiddle_c1 * i3) + (twiddle_s1 * r4);
    t2 = (twiddle_c1 * -r4) + (twiddle_s1 * i3);

    *r1_ptr++ = r3 + t1;
    *i1_ptr++ = t2 + i4;

    *r2_ptr-- = r3 - t1;
    *i2_ptr-- = t2 - i4;
  }
}

void do_small_real_fft(FFT_Split *input, t_uint fft_log2, long ifft)
{
  t_sample r1, r2, r3, r4, i1, i2;
  t_sample scale = 2.;

  t_sample *realp = input->realp;
  t_sample *imagp = input->imagp;

  if(fft_log2 < 1)
    return;

  if(ifft)
    scale = 1.;

  switch(fft_log2)
    {
    case 1:
      r1 = realp[0];
      r2 = imagp[0];

      realp[0] = (r1 + r2) * scale;
      imagp[0] = (r1 - r2) * scale;

      break;

    case 2:
      if(!ifft) {
	r3 = realp[0];
	r4 = realp[1];
	i1 = imagp[0];
	i2 = imagp[1];

	r1 = r3 + r4;
	r2 = r3 - r4;
	r3 = i1 + i2;
	r4 = i1 - i2;

	realp[0] = (r1 + r3) * 2;
	realp[1] = r2 * 2;
	imagp[0] = (r1 - r3) * 2;
	imagp[1] = -r4 * 2;
      }
      else {
	i1 = realp[0];
	r2 = realp[1] + realp[1];
	i2 = imagp[0];
	r4 = -imagp[1] - imagp[1];

	r1 = i1 + i2;
	r3 = i1 - i2;

	realp[0] = r1 + r2;
	realp[1] = r1 - r2;
	imagp[0] = r3 + r4;
	imagp[1] = r3 - r4;
      }
      break;
    }
}

void pass_trig_table(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass)
{
  t_uint size = 2 << pass;
  t_uint incr = size >> 1;
  t_uint loop = size;

  t_uint i;

  double r1, r2, i1, i2, r_0, i_0;
  double twiddle_c, twiddle_s;

  t_sample *r1_ptr = (t_sample *)input->realp;
  t_sample *i1_ptr = (t_sample *)input->imagp;
  t_sample *r2_ptr = (t_sample *)r1_ptr + incr;
  t_sample *i2_ptr = (t_sample *)i1_ptr + incr;

  for(i = 0; i < length; loop += size) {
    t_sample *tr_ptr = setup->tables[pass - PASS_TRIG_OFFSET].realp;
    t_sample *ti_ptr = setup->tables[pass - PASS_TRIG_OFFSET].imagp;

    for(; i < loop; i += 2) {
      twiddle_c = (double)*tr_ptr++;
      twiddle_s = (double)*ti_ptr++;

      r1 = *r1_ptr;
      i1 = *i1_ptr;
      r2 = *r2_ptr;
      i2 = *i2_ptr;

      r_0 = (r2 * twiddle_c) - (i2 * twiddle_s);
      i_0 = (r2 * twiddle_s) + (i2 * twiddle_c);

      *r1_ptr++ = (t_sample) (r1 + r_0);
      *i1_ptr++ = (t_sample) (i1 + i_0);

      *r2_ptr++ = (t_sample) (r1 - r_0);
      *i2_ptr++ = (t_sample) (i1 - i_0);
    }
    
    r1_ptr += incr;
    r2_ptr += incr;
    i1_ptr += incr;
    i2_ptr += incr;
  }
}

void pass_trig_table_reorder(FFT_Split *input, FFT_Setup *setup, t_uint length, t_uint pass)
{
  t_uint size = 2 << pass;
  t_uint incr = size >> 1;
  t_uint loop = size;
  t_uint offset = length >> (pass + 1);
  t_uint outer_loop = ((length >> 1) / size ) / ((t_uint) 1 << pass);
  t_uint i, j;
    
  t_sample r1, r2, r3, r4, r5;
  t_sample i1, i2, i3, i4, i5;
  t_sample twiddle_c, twiddle_s;

  t_sample *r1_ptr = input->realp;
  t_sample *i1_ptr = input->imagp;
  t_sample *r2_ptr = r1_ptr + incr;
  t_sample *i2_ptr = i1_ptr + incr;
  t_sample *r3_ptr = r1_ptr + offset;
  t_sample *i3_ptr = i1_ptr + offset;
  t_sample *r4_ptr = r3_ptr + incr;
  t_sample *i4_ptr = i3_ptr + incr;

  for(j = 0, i = 0; i < (length >> 1); loop += size) {
    t_sample *t_ptr_r = setup->tables[pass - PASS_TRIG_OFFSET].realp;
    t_sample *t_ptr_i = setup->tables[pass - PASS_TRIG_OFFSET].imagp;

    for(; i < loop; i += 2) {
      r1 = *r1_ptr;
      i1 = *i1_ptr;
      r2 = *r3_ptr;
      i2 = *i3_ptr;

      twiddle_c = *t_ptr_r++;
      twiddle_s = *t_ptr_i++;

      r5 = (r2 * twiddle_c) - (i2 * twiddle_s);
      i5 = (r2 * twiddle_s) + (i2 * twiddle_c);

      r3 = *r2_ptr;
      i3 = *i2_ptr;
      r4 = *r4_ptr;
      i4 = *i4_ptr;

      *r1_ptr++ = r1 + r5;
      *i1_ptr++ = i1 + i5;

      *r2_ptr++ = r1 - r5;
      *i2_ptr++ = i1 - i5;

      r5 = (r4 * twiddle_c) - (i4 * twiddle_s);
      i5 = (r4 * twiddle_s) + (i4 * twiddle_c);

      *r3_ptr++ = r3 + r5;
      *i3_ptr++ = i3 + i5;

      *r4_ptr++ = r3 - r5;
      *i4_ptr++ = i3 - i5;
    }

    r1_ptr += incr;
    r2_ptr += incr;
    r3_ptr += incr;
    r4_ptr += incr;
    i1_ptr += incr;
    i2_ptr += incr;
    i3_ptr += incr;
    i4_ptr += incr;

    if(!(++j % outer_loop)) {
      r1_ptr += offset;
      r2_ptr += offset;
      r3_ptr += offset;
      r4_ptr += offset;
      i1_ptr += offset;
      i2_ptr += offset;
      i3_ptr += offset;
      i4_ptr += offset;
    }
  }
}

void pass_3_reorder(FFT_Split *input, t_uint length, t_uint fft_log2)
{
  t_uint offset = (t_uint) 1 << (fft_log2 - 3);
  t_uint i, j, loop;

  t_sample sqrt_2_2 = (t_sample) SQRT_2_2;
  t_sample r1, r2, r3, r4, r5, r6, r7, r8;
  t_sample r9, r10, r11, r12, r13, r14, r15, r16;
  t_sample i1, i2, i3, i4, i5, i6, i7, i8;
  t_sample i9, i10, i11, i12, i13, i14, i15, i16;
  t_sample t1, t2, t3, t4;

  t_sample *r1_ptr = input->realp;
  t_sample *i1_ptr = input->imagp;

  t_sample *r2_ptr = r1_ptr + offset;
  t_sample *i2_ptr = i1_ptr + offset;

  for(j = 0, loop = length >> 6; j < 4; j++) {
    for(i=0; i<loop; i++) {
      r1  = *(r1_ptr + 0);
      r2  = *(r1_ptr + 1);
      r3  = *(r1_ptr + 2);
      r4  = *(r1_ptr + 3);

      r9  = *(r1_ptr + 4);
      r10 = *(r1_ptr + 5);
      r11 = *(r1_ptr + 6);
      r12 = *(r1_ptr + 7);

      r5  = *(r2_ptr + 0);
      r6  = *(r2_ptr + 1);
      r7  = *(r2_ptr + 2);
      r8  = *(r2_ptr + 3);

      r13 = *(r2_ptr + 4);
      r14 = *(r2_ptr + 5);
      r15 = *(r2_ptr + 6);
      r16 = *(r2_ptr + 7);

      i1  = *(i1_ptr + 0);
      i2  = *(i1_ptr + 1);
      i3  = *(i1_ptr + 2);
      i4  = *(i1_ptr + 3);

      i9  = *(i1_ptr + 4);
      i10 = *(i1_ptr + 5);
      i11 = *(i1_ptr + 6);
      i12 = *(i1_ptr + 7);

      i5  = *(i2_ptr + 0);
      i6  = *(i2_ptr + 1);
      i7  = *(i2_ptr + 2);
      i8  = *(i2_ptr + 3);

      i13 = *(i2_ptr + 4);
      i14 = *(i2_ptr + 5);
      i15 = *(i2_ptr + 6);
      i16 = *(i2_ptr + 7);

      t1 = sqrt_2_2 * (r6 + i6);
      t2 = sqrt_2_2 * (i8 - r8);
      t3 = sqrt_2_2 * (r14 + i14);
      t4 = sqrt_2_2 * (i16 - r16);

      *r1_ptr++ = r1 + r5;
      *r1_ptr++ = r2 + t1;
      *r1_ptr++ = r3 + i7;
      *r1_ptr++ = r4 + t2;
      *r1_ptr++ = r1 - r5;
      *r1_ptr++ = r2 - t1;
      *r1_ptr++ = r3 - i7;
      *r1_ptr++ = r4 - t2;

      *r2_ptr++ = r9  + r13;
      *r2_ptr++ = r10 + t3;
      *r2_ptr++ = r11 + i15;
      *r2_ptr++ = r12 + t4;
      *r2_ptr++ = r9  - r13;
      *r2_ptr++ = r10 - t3;
      *r2_ptr++ = r11 - i15;
      *r2_ptr++ = r12 - t4;

      t1 = sqrt_2_2 * (i6 - r6);
      t2 = -sqrt_2_2 * (r8 + i8);
      t3 = sqrt_2_2 * (i14 - r14);
      t4 = -sqrt_2_2 * (r16 + i16);

      *i1_ptr++ = i1 + i5;
      *i1_ptr++ = i2 + t1;
      *i1_ptr++ = i3 - r7;
      *i1_ptr++ = i4 + t2;
      *i1_ptr++ = i1 - i5;
      *i1_ptr++ = i2 - t1;
      *i1_ptr++ = i3 + r7;
      *i1_ptr++ = i4 - t2;

      *i2_ptr++ = i9  + i13;
      *i2_ptr++ = i10 + t3;
      *i2_ptr++ = i11 - r15;
      *i2_ptr++ = i12 + t4;
      *i2_ptr++ = i9  - i13;
      *i2_ptr++ = i10 - t3;
      *i2_ptr++ = i11 + r15;
      *i2_ptr++ = i12 - t4;
    }
    r1_ptr += offset;
    i1_ptr += offset;
    r2_ptr += offset;
    i2_ptr += offset;
  }
}

void pass_1_2_reorder(FFT_Split *input, t_uint length, t_uint fft_log2)
{
  t_uint i, j, loop;
  t_uint offset = length >> 1;

  t_sample *r1_ptr = input->realp;
  t_sample *i1_ptr = input->imagp;
  t_sample *r2_ptr = r1_ptr + offset;
  t_sample *i2_ptr = i1_ptr + offset;

  t_sample r1, r2, r3, r4, r5, r6, r7, r8;
  t_sample i1, i2, i3, i4, i5, i6, i7, i8;

  for(i=0; i< (length >> 2); i++) {
    r1 = *r1_ptr;
    r2 = *r2_ptr;
    r3 = *(r1_ptr + 1);
    r4 = *(r2_ptr + 1);

    *r1_ptr++ = r1 + r2;
    *r1_ptr++ = r1 - r2;
    *r2_ptr++ = r3 + r4;
    *r2_ptr++ = r3 - r4;

    i1 = *i1_ptr;
    i2 = *i2_ptr;
    i3 = *(i1_ptr + 1);
    i4 = *(i2_ptr + 1);

    *i1_ptr++ = i1 + i2;
    *i1_ptr++ = i1 - i2;
    *i2_ptr++ = i3 + i4;
    *i2_ptr++ = i3 - i4;
  }

  offset >>= 1;

  r1_ptr = input->realp;
  i1_ptr = input->imagp;
  r2_ptr = r1_ptr + offset;
  i2_ptr = i1_ptr + offset;

  for(j=0, loop = (length >> 4); j < 2; j++) {
    for(i=0;i<loop;i++) {
      r1 = *(r1_ptr + 0);
      r2 = *(r1_ptr + 1);
      r3 = *(r2_ptr + 0);
      r4 = *(r2_ptr + 1);
      r5 = *(r1_ptr + 2);
      r6 = *(r1_ptr + 3);
      r7 = *(r2_ptr + 2);
      r8 = *(r2_ptr + 3);

      i1 = *(i1_ptr + 0);
      i2 = *(i1_ptr + 1);
      i3 = *(i2_ptr + 0);
      i4 = *(i2_ptr + 1);
      i5 = *(i1_ptr + 2);
      i6 = *(i1_ptr + 3);
      i7 = *(i2_ptr + 2);
      i8 = *(i2_ptr + 3);

      *r1_ptr++ = r1 + r3;
      *r1_ptr++ = r2 + i4;
      *r1_ptr++ = r1 - r3;
      *r1_ptr++ = r2 - i4;
      *r2_ptr++ = r5 + r7;
      *r2_ptr++ = r6 + i8;
      *r2_ptr++ = r5 - r7;
      *r2_ptr++ = r6 - i8;

      *i1_ptr++ = i1 + i3;
      *i1_ptr++ = i2 - r4;
      *i1_ptr++ = i1 - i3;
      *i1_ptr++ = i2 + r4;
      *i2_ptr++ = i5 + i7;
      *i2_ptr++ = i6 - r8;
      *i2_ptr++ = i5 - i7;
      *i2_ptr++ = i6 + r8;      
    }

    r1_ptr += offset;
    i1_ptr += offset;
    r2_ptr += offset;
    i2_ptr += offset;
  }
}

void pass_3(FFT_Split *input, t_uint length, t_uint fft_log2)
{
  t_uint i;

  t_sample sqrt_2_2 = (t_sample) SQRT_2_2;
  t_sample r1, r2, r3, r4, r5, r6, r7, r8;
  t_sample i1, i2, i3, i4, i5, i6, i7, i8;
  t_sample t1, t2;

  t_sample *r1_ptr = input->realp;
  t_sample *i1_ptr = input->imagp;

  for(i=0; i < (length >> 3); i++) {
    r1 = *(r1_ptr + 0);
    r2 = *(r1_ptr + 1);
    r3 = *(r1_ptr + 2);
    r4 = *(r1_ptr + 3);
    r5 = *(r1_ptr + 4);
    r6 = *(r1_ptr + 5);
    r7 = *(r1_ptr + 6);
    r8 = *(r1_ptr + 7);

    i1 = *(i1_ptr + 0);
    i2 = *(i1_ptr + 1);
    i3 = *(i1_ptr + 2);
    i4 = *(i1_ptr + 3);
    i5 = *(i1_ptr + 4);
    i6 = *(i1_ptr + 5);
    i7 = *(i1_ptr + 6);
    i8 = *(i1_ptr + 7);

    t1 = sqrt_2_2 * (r6 + i6);
    t2 = sqrt_2_2 * (i8 - r8);

    *r1_ptr++ = r1 + r5;
    *r1_ptr++ = r2 + t1;
    *r1_ptr++ = r3 + i7;
    *r1_ptr++ = r4 + t2;
    *r1_ptr++ = r1 - r5;
    *r1_ptr++ = r2 - t1;
    *r1_ptr++ = r3 - i7;
    *r1_ptr++ = r4 - t2;

    t1 = sqrt_2_2 * (i6 - r6);
    t2 = -sqrt_2_2 * (r8 + i8);

    *i1_ptr++ = i1 + i5;
    *i1_ptr++ = i2 + t1;
    *i1_ptr++ = i3 - r7;
    *i1_ptr++ = i4 + t2;
    *i1_ptr++ = i1 - i5;
    *i1_ptr++ = i2 - t1;
    *i1_ptr++ = i3 + r7;
    *i1_ptr++ = i4 - t2;
  }
}

void do_small_fft(FFT_Split *input, t_uint fft_log2)
{
  t_sample r1, r2, r3, r4, r5, r6, r7, r8;
  t_sample i1, i2, i3, i4, i5, i6, i7, i8;
  t_sample t1, t2, t3, t4;

  t_sample *realp = input->realp;
  t_sample *imagp = input->imagp;

  if(fft_log2 < 1)
    return;

  switch(fft_log2)
    {
    case 1:
      
      r1 = realp[0];
      r2 = realp[1];
      i1 = imagp[0];
      i2 = imagp[1];

      realp[0] = r1 + r2;
      realp[1] = r1 - r2;
      imagp[0] = i1 + i2;
      imagp[1] = i1 - i2;

      break;

    case 2:

      r5 = realp[0];
      r6 = realp[1];
      r2 = realp[2];
      r4 = realp[3];
      i5 = imagp[0];
      i6 = imagp[1];
      i2 = imagp[2];
      i4 = imagp[3];

      r1 = r5 + r2;
      r2 = r5 - r2;
      r3 = r6 + r4;
      r4 = r6 - r4;
      i1 = i5 + i2;
      i2 = i5 - i2;
      i3 = i6 + i4;
      i4 = i6 - i4;

      realp[0] = r1 + r3;
      realp[1] = r2 + i4;
      realp[2] = r1 - r3;
      realp[3] = r2 - i4;
      imagp[0] = i1 + i3;
      imagp[1] = i2 - r4;
      imagp[2] = i1 - i3;
      imagp[3] = i2 + i4;

      break;

    case 3:

      t1 = realp[0];
      t3 = realp[1];
      t2 = realp[2];
      t4 = realp[3];
      r2 = realp[4];
      r6 = realp[5];
      r4 = realp[6];
      r8 = realp[7];

      r1 = t1 + r2;
      r2 = t1 - r2;
      r3 = t2 + r4;
      r4 = t2 - r4;
      r5 = t3 + r6;
      r6 = t3 - r6;
      r7 = t4 + r8;
      r8 = t4 - r8;

      t1 = imagp[0];
      t3 = imagp[1];
      t2 = imagp[2];
      t4 = imagp[3];
      i2 = imagp[4];
      i6 = imagp[5];
      i4 = imagp[6];
      i8 = imagp[7];

      i1 = t1 + i2;
      i2 = t1 - i2;
      i3 = t2 + i4;
      i4 = t2 - i4;
      i5 = t3 + i6;
      i6 = t3 - i6;
      i7 = t4 + i8;
      i8 = t4 - i8;

      realp[0] = r1 + r3;
      realp[1] = r2 + i4;
      realp[2] = r1 - r3;
      realp[3] = r2 - i4;
      realp[4] = r5 + r7;
      realp[5] = r6 + i8;
      realp[6] = r5 - r7;
      realp[7] = r6 - i8;

      imagp[0] = i1 + i3;
      imagp[1] = i2 - r4;
      imagp[2] = i1 - i3;
      imagp[3] = i2 + r4;
      imagp[4] = i5 + i7;
      imagp[5] = i6 - r8;
      imagp[6] = i5 - i7;
      imagp[7] = i6 + r8;

      pass_3(input,(t_uint)8,(t_uint)3);

      break;
    }
}

t_uint int_log2(t_uint in, t_uint *inexact)
{
  t_uint temp = in;
  t_uint out = 0;

  if(!in)
    return 0;

  while(temp) {
    temp >>= 1;
    out++;
  }

  if(in == (t_uint) 1 << (out - (t_uint)1)) {
    out--;
    *inexact = 0;
  }
  else
    *inexact = 1;

  return out;
}

t_uint calculate_fft_size(t_uint input_size, t_uint *fft_size_log2)
{
  t_uint inexact = 0;
  t_uint calc_temp = int_log2(input_size, &inexact);

  *fft_size_log2 = calc_temp;

  if(!input_size)
    return 0;

  return (t_uint) 1 << calc_temp;
}

void fft_fill_table32(FFT_Split32 *table, t_uint length)
{
  t_uint i;
  t_float32 *table_cos = table->realp;
  t_float32 *table_sin = table->imagp;

  for(i=0; i<length; i++) {
    double angle = -((double) i) * M_PI / (double) length;

    *table_cos++ = (t_float32)cos(angle);
    *table_sin++ = (t_float32)sin(angle);
  }
}

FFT_Setup32 *create_setup32(t_uint max_fft_log2)
{
  FFT_Setup32 *setup = (FFT_Setup32 *)aligned_getbytes(sizeof(FFT_Setup32));
  t_uint i;

  SSE_Exists = SSE2_check();
  
  for(i=FFTLOG2_TRIG_OFFSET; i<=max_fft_log2; i++) {
    t_uint length = (t_uint) 1 << (i-1);

    setup->tables[i - FFTLOG2_TRIG_OFFSET].realp =
      (t_float32 *)aligned_getbytes(sizeof(t_float32) * 2 * length);
    setup->tables[i - FFTLOG2_TRIG_OFFSET].imagp =
      setup->tables[i - FFTLOG2_TRIG_OFFSET].realp + length;
    fft_fill_table32(&setup->tables[i - FFTLOG2_TRIG_OFFSET], length);
  }

  setup->max_fft_log2 = max_fft_log2;

  return setup;
}

void destroy_setup32(FFT_Setup32 *setup)
{ 
  t_uint max_fft_log2 = setup->max_fft_log2;
  t_uint i;

  for(i=FFTLOG2_TRIG_OFFSET; i<=max_fft_log2; i++) {
    t_uint length = (t_uint) 1 << (i-1);
    aligned_freebytes(setup->tables[i-FFTLOG2_TRIG_OFFSET].realp,
	      sizeof(t_float32)*2*length);
  }
  aligned_freebytes(setup, sizeof(FFT_Setup32));
}

void fft_fill_table(FFT_Split *table, t_uint length)
{
  t_uint i;
  t_sample *table_cos = table->realp;
  t_sample *table_sin = table->imagp;

  for(i=0; i<length; i++) {
    double angle = -((double) i) * M_PI / (double) length;

    *table_cos++ = (t_sample)cos(angle);
    *table_sin++ = (t_sample)sin(angle);
  }
}

FFT_Setup *create_setup(t_uint max_fft_log2)
{
  FFT_Setup *setup = (FFT_Setup *)aligned_getbytes(sizeof(FFT_Setup));
  t_uint i;

  SSE_Exists = SSE2_check();
  
  for(i=FFTLOG2_TRIG_OFFSET; i<=max_fft_log2; i++) {
    t_uint length = (t_uint) 1 << (i-1);

    setup->tables[i - FFTLOG2_TRIG_OFFSET].realp =
      (t_sample *)aligned_getbytes(sizeof(t_sample) * 2 * length);
    setup->tables[i - FFTLOG2_TRIG_OFFSET].imagp =
      setup->tables[i - FFTLOG2_TRIG_OFFSET].realp + length;
    fft_fill_table(&setup->tables[i - FFTLOG2_TRIG_OFFSET], length);
  }

  setup->max_fft_log2 = max_fft_log2;

  return setup;
}

void destroy_setup(FFT_Setup *setup)
{ 
  t_uint max_fft_log2 = setup->max_fft_log2;
  t_uint i;

  for(i=FFTLOG2_TRIG_OFFSET;i<=max_fft_log2;i++) {
    t_uint length = (t_uint) 1 << (i-1);
    aligned_freebytes(setup->tables[i-FFTLOG2_TRIG_OFFSET].realp,
	      sizeof(t_sample)*2*length);
  }
  aligned_freebytes(setup, sizeof(FFT_Setup));
}

void fft_internal32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2)
{
  t_uint_fft length = (t_uint_fft) 1 << fft_log2;
  t_uint_fft i;
  
  if(fft_log2 < 4) {
    do_small_fft32(input, fft_log2);
    return;
  }
#ifdef VECTOR_F64_128BIT
  if((t_uint_fft) input->realp % 16 || (t_uint_fft) input->imagp % 16 || !SSE_Exists)
#endif
    {
  
      pass_1_2_reorder32(input, length, fft_log2);

      if(fft_log2 > 5) 
      	pass_3_reorder32(input, length, fft_log2);
      else
      	pass_3_32(input, length, fft_log2);
  
      for(i=3; i < (fft_log2 >> 1); i++)
      	pass_trig_table_reorder32(input, setup, length, i);
      
      for(; i < fft_log2; i++)
      	pass_trig_table32(input, setup, length, i);
    }
#ifdef VECTOR_F64_128BIT
  else
    {
      
      
      pass_1_2_reorder_simd32(input, length);

      if(fft_log2 > 5)
	pass_3_reorder_simd32(input, length);
      else
	pass_3_32(input, length, fft_log2);

      for(i=3; i<(fft_log2 >> 1); i++)
	pass_trig_table_reorder_simd32(input, setup, length, i);

      for(; i<fft_log2;i++)
	pass_trig_table_simd32(input, setup, length, i);
    }
#endif
}

void fft_internal(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2)
{
  t_uint_fft length = (t_uint_fft) 1 << fft_log2;
  t_uint_fft i;
  
  if(fft_log2 < 4) {
    do_small_fft(input, fft_log2);
    return;
  }
#ifdef VECTOR_F64_128BIT
  if((t_uint_fft) input->realp % 16 || (t_uint_fft) input->imagp % 16 || !SSE_Exists)
#endif
    {
  
      pass_1_2_reorder(input, length, fft_log2);

      if(fft_log2 > 5) 
      	pass_3_reorder(input, length, fft_log2);
      else
      	pass_3(input, length, fft_log2);
  
      for(i=3; i < (fft_log2 >> 1); i++)
      	pass_trig_table_reorder(input, setup, length, i);
      
      for(; i < fft_log2; i++)
      	pass_trig_table(input, setup, length, i);
    }
#ifdef VECTOR_F64_128BIT
  else
    {
      
      
      pass_1_2_reorder_simd(input, length);

      if(fft_log2 > 5)
	pass_3_reorder_simd(input, length);
      else
	pass_3(input, length, fft_log2);

      for(i=3; i<(fft_log2 >> 1); i++)
	pass_trig_table_reorder_simd(input, setup, length, i);

      for(; i<fft_log2;i++)
	pass_trig_table_simd(input, setup, length, i);
    }
#endif
}

void do_fft32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2)
{
  fft_internal32(input, setup, fft_log2);
}

void do_ifft32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2)
{
  FFT_Split32 swap;

  swap.realp = input->imagp;
  swap.imagp = input->realp;

  fft_internal32(&swap, setup, fft_log2);
}

void do_real_fft32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2)
{
  if(fft_log2 < 3) {
    do_small_real_fft32(input, fft_log2, 0L);
    return;
  }

  do_fft32(input, setup, fft_log2 - 1);
  pass_real_trig_table32(input, setup, fft_log2, 0L);
}

void do_real_ifft32(FFT_Split32 *input, FFT_Setup32 *setup, t_uint fft_log2)
{
  if(fft_log2 < 3) {
    do_small_real_fft32(input, fft_log2, 1L);
    return;
  }

  pass_real_trig_table32(input, setup, fft_log2, 1L);
  do_ifft32(input, setup, fft_log2 - 1);
}

void do_fft(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2)
{
  fft_internal(input, setup, fft_log2);
}

void do_ifft(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2)
{
  FFT_Split swap;

  swap.realp = input->imagp;
  swap.imagp = input->realp;

  fft_internal(&swap, setup, fft_log2);
}

void do_real_fft(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2)
{
  if(fft_log2 < 3) {
    do_small_real_fft(input, fft_log2, 0L);
    return;
  }

  do_fft(input, setup, fft_log2 - 1);
  pass_real_trig_table(input, setup, fft_log2, 0L);
}

void do_real_ifft(FFT_Split *input, FFT_Setup *setup, t_uint fft_log2)
{
  if(fft_log2 < 3) {
    do_small_real_fft(input, fft_log2, 1L);
    return;
  }

  pass_real_trig_table(input, setup, fft_log2, 1L);
  do_ifft(input, setup, fft_log2 - 1);
}

void unzip_complex32(t_float32 *input, FFT_Split32 *output, t_uint half_length)
{
  t_uint i;
  
  t_float32 *realp = output->realp;
  t_float32 *imagp = output->imagp;

  for(i=0; i<half_length; i++) {
    *realp++ = *input++;
    *imagp++ = *input++;
  }
}

void unzip_complex(t_sample *input, FFT_Split *output, t_uint half_length)
{
  t_uint i;
  
  t_sample *realp = output->realp;
  t_sample *imagp = output->imagp;

  for(i=0; i<half_length; i++) {
    *realp++ = *input++;
    *imagp++ = *input++;
  }
}

void unzip_zero(t_sample *input, FFT_Split *output, t_uint in_length, t_uint log2n)
{
  t_uint i;
  t_sample temp = 0.;

  t_sample *realp = output->realp;
  t_sample *imagp = output->imagp;

  t_uint half_length;
  
  if(((t_uint) 1 << log2n) < in_length) {
    in_length = (t_uint) 1 << log2n;
  }
  if(in_length & 1) {
    temp = input[in_length - 1];
  }

  half_length = in_length >> (t_uint)1;

  unzip_complex(input, output, half_length);

  if(((t_uint) 1 << log2n) > in_length) {
    
    realp[half_length] = temp;
    imagp[half_length] = 0.;

    for(i = (in_length >> (t_uint) 1) + 1;
	i < ((t_uint) 1 << (log2n - (t_uint) 1));
	i++) {
      
      realp[i] = 0.;
      imagp[i] = 0.;
    }
  }
}

void zip_sample32(FFT_Split32 *input, t_float32 *output, t_uint half_length)
{
  t_uint i;

  t_float32 *realp = input->realp;
  t_float32 *imagp = input->imagp;

  for(i=0; i<half_length; i++) {
    *output++ = *realp++;
    *output++ = *imagp++;
  }
}

void zip_sample(FFT_Split *input, t_sample *output, t_uint half_length)
{
  t_uint i;

  t_sample *realp = input->realp;
  t_sample *imagp = input->imagp;

  for(i=0; i<half_length; i++) {
    *output++ = *realp++;
    *output++ = *imagp++;
  }
}
