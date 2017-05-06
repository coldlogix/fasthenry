#ifndef __joelself_H__
#define __joelself_H__

#include "induct.h"

/* joelself.c */
double self(double, double, double);
int edges_parallel(FILAMENT*, FILAMENT*, double*, int*);
void get_wid(FILAMENT*, double*);
void get_height(FILAMENT*, double*, double*);
double exact_mutual(FILAMENT*, FILAMENT*, int, double*, double*,
    enum degen_type, enum degen_type);
// void fill_4(double*, double, double, double);
// double eval_eq(double, double, double, double);
// double log_term(double, double, double, double, double);
// double tan_term(double, double, double, double, double);
int lookup(FILAMENT*, FILAMENT*, int, double*, double*, double*, double*,
    int*, Table***, int*);
// void fill_dims(FILAMENT*, FILAMENT*, double*, double*, double*, int);
// void fill_dims_seg(FILAMENT*, FILAMENT*, double*, double*, double*, int);
// int find_dims(double*, int, Table**, double*, int*, Table***);
void put_in_table(FILAMENT*, FILAMENT*, int, double, double*, int, Table**,
    int);
void init_table(void);
int get_table_mem(void);
void destroy_table(void);
// char *AllocAnEntry(AllocInfo*);
// void DestroyEntries(AllocInfo*);
// int MemoryForEntries(AllocInfo*);
// double brick_to_brick(double, double, double, double, double, double, double,
//     double, double);
// double flat_to_flat_tape(double, double, double, double, double, double,
//     double);
// double eval_eq_tape(double, double, double, double);
// double flat_to_skinny_tape(double, double, double, double, double, double,
//     double);
// double eval_eq_tape2(double, double, double, double);
// double tan_tape(double, double, double, double);
// double tape_to_fil(double, double, double, double, double, double);
// double brick_to_fil(double, double, double, double, double, double, double);

#endif
