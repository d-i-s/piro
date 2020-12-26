#include "m_types.h"

/*
 *  we need to redefine memory allocation to consider memory alignment for win
 *	
 *
 *  Copyright (c) 1997-1999 Miller Puckette. All rights reserved.
 *  Copyright (c) 2019 Marco Matteo Markidis. All rights reserved.
 */

void *aligned_getbytes(size_t nbytes)
{
  void *ret;
#if defined(UNIX) || defined (MACOSX)
  ret = getbytes(nbytes);
#elif defined (NT)

  if (nbytes < 1) nbytes = 1;
  ret = (void *)ALIGNED_MALLOC(nbytes);
  if(!ret)
    post("pd: getbytes() failed -- out of memory");
  memset(ret, 0, nbytes);
#endif
  return (ret);
}

void *aligned_resizebytes(void *old, size_t oldsize, size_t newsize)
{
  void *ret;
#if defined(UNIX) || defined (MACOSX)
  ret = resizebytes(old, oldsize, newsize);
#elif defined (NT)
  if (newsize < 1) newsize = 1;
  if (oldsize < 1) oldsize = 1;
  ret = (void *)ALIGNED_REALLOC((char *)old, newsize);
  if (newsize > oldsize && ret)
    memset(((char *)ret) + oldsize, 0, newsize - oldsize);
  if (!ret)
    post("pd: resizebytes() failed -- out of memory");
#endif
  return (ret);
}

void aligned_freebytes(void *fatso, size_t nbytes)
{
#if defined(UNIX) || defined (MACOSX)
  freebytes(fatso, nbytes);
#elif defined (NT)
  if (nbytes == 0)
    nbytes = 1;
  ALIGNED_FREE(fatso);
#endif
}
