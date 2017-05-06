#ifndef __barnoldi_H__
#define __barnoldi_H__

#include "induct.h"

typedef int (*barnoldi_cb)(double**, ssystem*, double**, charge*, SYS*, int,
    int, int);
int ArnoldiROM(double**, double**, double**, char**, int, int, int, int,
    barnoldi_cb, SYS*, ssystem*, charge*);
// int qr(double**, double**, double**, int, int, int);
// int qr_P(double**, double**, double**, double**, int, int, int, char*);
int dumpROM(FILE*, double**, double**, double**, double**, int, int, int);
void dumpROMequiv_circuit(FILE*, double**, double**, double**, double**,
    int, int, int, char*, char*, SYS*);
// int dumpROMbin(FILE*, double**, double**, double**, double**, int, int, int);
int createMRMt(char**, SYS*);
// int createMRMtinvMLMt(double***, SYS*, char*);
int realComputePsi(double**, ssystem*, double**, charge*, SYS*, int, int, int);
int realMatVect(double**, ssystem*, double**, charge*, SYS*, int, int, int);
// int printRowCol(double**, int, int, int, int);
void formMLMt(SYS*);
// void ZeroMatrix(double**, int, int);

#endif
