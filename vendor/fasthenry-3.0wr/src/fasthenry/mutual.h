#ifndef __mutal_H__
#define __mutal_H__

#include "induct.h"

/* mutual.c */
double mutual(FILAMENT*, FILAMENT*);
// void print_infinity_warning(FILAMENT*, FILAMENT*);
// void findfourfils(FILAMENT*, FILAMENT*);
double selfterm(FILAMENT*);
double mutualfil(FILAMENT*, FILAMENT*);
// double magdiff2(FILAMENT*, int, FILAMENT*, int);
// double mut_rect(double, double);
double dotprod(FILAMENT*, FILAMENT*);
double fourfil(FILAMENT*, FILAMENT*);
// double parallel_fils(FILAMENT*, FILAMENT*, int, double*, double*, double);

#endif

