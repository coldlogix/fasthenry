#ifndef __readGeom_H__
#define __readGeom_H__

#include <stdio.h>
#include "induct.h"

/* readGeom.c */
int readGeom(FILE*, SYS*);
// int dodot(char*, SYS*);
// int changeunits(char*, SYS*);
// int addexternal(char*, SYS*);
// int choosefreqs(char*, SYS*);
// int old_equivnodes(char*, SYS*);
// int dodefault(char*);
// int addnode(char*, SYS*, NODES**, int);
NODES *makenode(char*, int, double, double, double, int, SYS*, int);
// int addseg(char*, SYS*, int, SEGMENT**);
SEGMENT *makeseg(char*, NODES*, NODES*, double, double, double,
#if SUPERCON == ON
    double,
#endif
    int, int, double, double, double*, int, int, SYS*, int);
// int addgroundplane(char*, SYS*, GROUNDPLANE**);
// int nothing(char*);
// char *getaline(FILE*);
// char *plusline(FILE*);
// char *getoneline(FILE*);
// void savealine(char*);
int notblankline(char*);
void tolowercase(char*);
int is_nonuni_gp(GROUNDPLANE*);

#endif
