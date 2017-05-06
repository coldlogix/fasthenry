#ifndef __dist_betw_fils_H__
#define __dist_betw_fils_H__

#include "induct.h"

/* dist_betw_fils.c */
double dist_betw_fils(FILAMENT*, FILAMENT*, int*);
// void getD(FILAMENT*, double*);
// void getr(double*, double*, double*, double*, double, double*);
double vdotp(double*, double*);
// double dist_between(double, double, double, double, double, double);
// double min_endpt_sep(FILAMENT*, FILAMENT*);
// double dist_betw_pt_and_fil(FILAMENT*, double*, double*, double, FILAMENT*,
//     double);
// double aspectratio(FILAMENT*);
void fill_Gquad(void);
// void findnfils(FILAMENT*, FILAMENT*, int);
// void gquad_weights(int, double*, double*);


#endif
