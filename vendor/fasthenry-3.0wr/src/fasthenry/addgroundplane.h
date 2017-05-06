#ifndef __addgroundplane_H__
#define __addgroundplane_H__

#include "induct.h"

int checkmiddlepoint(double*, double*, double*, int, int, int);
int checkplaneformula(double*, double*, double*, double, double, double,
    int, int, int);
double findsegmentwidth(double*, double*, double*, int, int, int, int);
void doincrement(double, double, double, double, double, double, int,
    double*, double*, double*);
void dounitvector(double, double, double, double, double, double,
    double*, double*, double*);
void make_nodelist(NODELIST*, char*, double, double, double);
SPATH *path_through_gp(SYS*, NODES*, NODES*, GROUNDPLANE*);
void clear_marks(SYS*);
void dump_mesh_coords(SYS*);
void dump_ascii_mesh_coords(SYS*);

#endif
