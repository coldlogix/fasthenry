#ifndef __cx_ludecomp_H__
#define __cx_ludecomp_H__

#include "cmplx.h"

/* cx_ludecomp.c */
CX **cx_ludecomp(CX**, int, int);
void cx_lu_solve(CX**, CX*, CX*, int);

#endif
