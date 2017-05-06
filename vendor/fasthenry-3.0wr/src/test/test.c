/* This is the main part of the code */

#include "../fasthenry/sparse/spMatrix.h"
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <ctype.h>
#include "../fasthenry/cmplx.h"

#ifdef MPI
#include <mpi.h>
#endif


int main(int argc, char **argv)
{
  char* spMatrix;
  FILE* fop;
  CX *Ima, *Imb;
  unsigned long size;

  spMatrix=NULL;
  Ima= NULL;

  /* Load spMatrix */
  fop=fopen("spMatrix.bin","rt");
  if (fop==NULL)
  {
    fprintf(stderr, "Can not open spMatrix.bin file\n");
    return -1;
  }
  spMatrixLoad(&spMatrix, fop);
  fclose (fop);

  /* Load vector */
  fop=fopen("Ima_0.bin","rb");
  if (fop==NULL)
  {
    fprintf(stderr, "Can not open Ima_0.bin file\n");
    return -1;
  }
  {
    fpos_t position;
    fseek(fop,0, SEEK_END);
    fgetpos(fop, &position);
    fseek(fop,0, SEEK_SET);

    size=(unsigned long)position.__pos;
    Ima=(CX*)malloc(size);
    fread(Ima, size, 1,fop);

  }
  fclose (fop);

  /* Load vector */
  fop=fopen("Imb_0.bin","rb");
  if (fop==NULL)
  {
    fprintf(stderr, "Can not open Imb_0.bin file\n");
    return -1;
  }
  {
    Imb=(CX*)malloc(size);
    fread(Imb, size, 1,fop);
  }
  fclose (fop);

    printf("Matrix size %u\n", size/sizeof(CX));

  /* Do the magic */
  spSolve(spMatrix, (spREAL*)Ima, (spREAL*)Ima);

  /* Check to result to reference */
  {
    unsigned int i,error;
    char* dataa, *datab;
    dataa=(char*)Ima;
    datab=(char*)Imb;
    error=0;
    for (i=0; (i<size) && (error==0); i++)
    {
       if ((*dataa)!=(*datab))
       {
         error=1;
       }
       else
       {
         dataa++;
         datab++;
       }
    }

    if (error==0)
    {
      printf("Test successful\n");
    }
    else
    {
      printf("Test failed\n");
      return -1;
    }
  }
  /* Clean up */
  free (Ima); Ima=NULL;
  free (Imb); Imb=NULL;
  spDestroy(spMatrix);

}
