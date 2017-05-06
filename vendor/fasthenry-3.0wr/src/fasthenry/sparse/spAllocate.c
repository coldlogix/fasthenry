/*
 *  MATRIX ALLOCATION MODULE
 *
 *  Author:                     Advising professor:
 *      Kenneth S. Kundert          Alberto Sangiovanni-Vincentelli
 *      UC Berkeley
 *
 *  This file contains the allocation and deallocation routines for the
 *  sparse matrix routines.
 *
 *  >>> User accessible functions contained in this file:
 *  spCreate
 *  spDestroy
 *  spError
 *  spWhereSingular
 *  spGetSize
 *  spSetReal
 *  spSetComplex
 *  spFillinCount
 *  spElementCount
 *
 *  >>> Other functions contained in this file:
 *  spcGetElement
 *  spcGetFillin
 */


/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and
 *  its documentation for any purpose and without fee is hereby granted,
 *  provided that the copyright notices appear in all copies and
 *  supporting documentation and that the authors and the University of
 *  California are properly credited.  The authors and the University of
 *  California make no representations as to the suitability of this
 *  software for any purpose.  It is provided `as is', without express
 *  or implied warranty.
 */

#ifndef lint
static char copyright[] =
    "Sparse1.3: Copyright (c) 1985,86,87,88 by Kenneth S. Kundert";
static char RCSid[] =
    "@(#)$Header: spAllocate.c,v 1.3 88/06/24 05:00:11 kundert Exp $";
#endif



/*
 *  IMPORTS
 *
 *  >>> Import descriptions:
 *  spConfig.h
 *      Macros that customize the sparse matrix routines.
 *  spMatrix.h
 *      Macros and declarations to be imported by the user.
 *  spDefs.h
 *      Matrix type and macro definitions for the sparse matrix routines.
 */

#define spINSIDE_SPARSE
#include "spConfig.h"
#include "spMatrix.h"
#include "spDefs.h"

#ifdef MPI
#include <mpi.h>
#endif


static void InitializeElementBlocks(MatrixPtr, int, int );
static void AllocateBlockOfAllocationList( MatrixPtr );

void spAllocateBlockOfAllocationHashList( MatrixPtr Matrix );
void spAllocateBlockOfAllocationList( MatrixPtr Matrix );
void spRecordAllocationFunc( MatrixPtr Matrix, char *AllocatedPtr, unsigned long AllocatedSize, unsigned int AllocatedType, const char* filename, unsigned int linenumber );
void spRemoveRecordAllocationFunc (MatrixPtr Matrix, char *AllocatedPtr, const char* filename, unsigned int linenumber);
spAllocationListHashPtr spFindRecordAllocation (MatrixPtr Matrix, char *Ptr );
unsigned long spGetMemoryBlockLength (MatrixPtr Matrix, char* AllocatedPtr);



/*
 *  MATRIX ALLOCATION
 *
 *  Allocates and initializes the data structures associated with a matrix.
 *
 *  >>> Returned:
 *  A pointer to the matrix is returned cast into the form of a pointer to
 *  a character.  This pointer is then passed and used by the other matrix
 *  routines to refer to a particular matrix.  If an error occurs, the NULL
 *  pointer is returned.
 *
 *  >>> Arguments:
 *  Size  <input>  (int)
 *      Size of matrix or estimate of size of matrix if matrix is EXPANDABLE.
 *  Complex  <input>  (int)
 *      Type of matrix.  If Complex is 0 then the matrix is real, otherwise
 *      the matrix will be complex.  Note that if the routines are not set up
 *      to handle the type of matrix requested, then a spPANIC error will occur.
 *      Further note that if a matrix will be both real and complex, it must
 *      be specified here as being complex.
 *  pError  <output>  (int *)
 *      Returns error flag, needed because function spError() will not work
 *      correctly if spCreate() returns NULL.
 *
 *  >>> Local variables:
 *  AllocatedSize  (int)
 *      The size of the matrix being allocated.
 *  Matrix  (MatrixPtr)
 *      A pointer to the matrix frame being created.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 *  spPANIC
 *  Error is cleared in this routine.
 */

char *
spCreate( int Size, BOOLEAN Complex, int *pError )

{
register  unsigned  SizePlusOne;
register  MatrixPtr  Matrix;
register  int  I;
int  AllocatedSize;

/* Begin `spCreate'. */
/* Clear error flag. */
    *pError = spOKAY;

/* Test for valid size. */
    if ((Size < 0) OR (Size == 0 AND NOT EXPANDABLE))
    {   *pError = spPANIC;
        return NULL;
    }

/* Test for valid type. */
#if NOT spCOMPLEX
    if (Complex)
    {   *pError = spPANIC;
        return NULL;
    }
#endif
#if NOT REAL
    if (NOT Complex)
    {   *pError = spPANIC;
        return NULL;
    }
#endif

/* Create Matrix. */
    AllocatedSize = MAX( Size, MINIMUM_ALLOCATED_SIZE );
    SizePlusOne = (unsigned)(AllocatedSize + 1);

    CALLOC(Matrix, struct MatrixFrame, 1);
    if (Matrix == NULL)
    {
        *pError = spNO_MEMORY;
        return NULL;
    }

/* Initialize matrix */
    Matrix->ID = SPARSE_ID;
    Matrix->Complex = Complex;
    Matrix->PreviousMatrixWasComplex = Complex;
    Matrix->Factored = NO;
    Matrix->Elements = 0;
    Matrix->Error = *pError;
    Matrix->Fillins = 0;
    Matrix->Reordered = NO;
    Matrix->NeedsOrdering = YES;
    Matrix->NumberOfInterchangesIsOdd = NO;
    Matrix->Partitioned = NO;
    Matrix->RowsLinked = NO;
    Matrix->InternalVectorsAllocated = NO;
    Matrix->SingularCol = 0;
    Matrix->SingularRow = 0;
    Matrix->Size = Size;
    Matrix->AllocatedSize = AllocatedSize;
    Matrix->ExtSize = Size;
    Matrix->AllocatedExtSize = AllocatedSize;
    Matrix->CurrentSize = 0;
    Matrix->ExtToIntColMap = NULL;
    Matrix->ExtToIntRowMap = NULL;
    Matrix->IntToExtColMap = NULL;
    Matrix->IntToExtRowMap = NULL;
    Matrix->MarkowitzRow = NULL;
    Matrix->MarkowitzCol = NULL;
    Matrix->MarkowitzProd = NULL;
    Matrix->DoCmplxDirect = NULL;
    Matrix->DoRealDirect = NULL;
    Matrix->Intermediate = NULL;
    Matrix->RelThreshold = DEFAULT_THRESHOLD;
    Matrix->AbsThreshold = 0.0;

    Matrix->TopOfAllocationList = NULL;
    Matrix->RecordsRemaining = 0;
    Matrix->ElementsRemaining = 0;
    Matrix->FillinsRemaining = 0;

    RecordAllocation( Matrix, (char *)Matrix, sizeof(struct MatrixFrame), AllocTypeMatrixFrame );
    if (Matrix->Error == spNO_MEMORY) goto MemoryError;

/* Take out the trash. */
    Matrix->TrashCan.Real = 0.0;
#if spCOMPLEX
    Matrix->TrashCan.Imag = 0.0;
#endif
    Matrix->TrashCan.Row = 0;
    Matrix->TrashCan.Col = 0;
    Matrix->TrashCan.NextInRow = NULL;
    Matrix->TrashCan.NextInCol = NULL;
#if INITIALIZE
    Matrix->TrashCan.pInitInfo = NULL;
#endif

/* Allocate space in memory for Diag pointer vector. */
    spAlloc( Matrix, &Matrix->Diag, sizeof(ElementPtr)*SizePlusOne, AllocTypePtrArray );
    if (Matrix->Diag == NULL)
        goto MemoryError;

/* Allocate space in memory for FirstInCol pointer vector. */
    spAlloc( Matrix, &Matrix->FirstInCol, sizeof(ElementPtr)*SizePlusOne, AllocTypePtrArray );
    if (Matrix->FirstInCol == NULL)
        goto MemoryError;

/* Allocate space in memory for FirstInRow pointer vector. */
    spAlloc( Matrix, &Matrix->FirstInRow, sizeof(ElementPtr)*SizePlusOne, AllocTypePtrArray );
    if (Matrix->FirstInRow == NULL)
        goto MemoryError;

/* Allocate space in memory for IntToExtColMap vector. */
    spAlloc( Matrix, &Matrix->IntToExtColMap, sizeof(int)*SizePlusOne, AllocTypeGeneric );
    if ( Matrix->IntToExtColMap == NULL)
        goto MemoryError;

/* Allocate space in memory for IntToExtRowMap vector. */
    spAlloc( Matrix, &Matrix->IntToExtRowMap, sizeof(int)*SizePlusOne, AllocTypeGeneric );
    if ( Matrix->IntToExtRowMap  == NULL)
        goto MemoryError;


/* Initialize MapIntToExt vectors. */
    for (I = 1; I <= AllocatedSize; I++)
    {
        Matrix->IntToExtRowMap[I] = I;
        Matrix->IntToExtColMap[I] = I;
    }

    /* SRW */
#if BUILDHASH
    Matrix->ElementHashTab = NULL;
#endif
#if BITFIELD
    Matrix->BitField = NULL;
#endif

#if TRANSLATE
    /* Allocate space in memory for ExtToIntColMap vector. */
    spAlloc( Matrix, &Matrix->ExtToIntColMap, sizeof(int)*SizePlusOne, AllocTypeGeneric );
    if ( Matrix->ExtToIntColMap  == NULL)
        goto MemoryError;


    /* Allocate space in memory for ExtToIntRowMap vector. */
    spAlloc( Matrix, &Matrix->ExtToIntRowMap, sizeof(int)*SizePlusOne, AllocTypeGeneric );
    if ( Matrix->ExtToIntRowMap  == NULL)
        goto MemoryError;


    /* Initialize MapExtToInt vectors. */
    for (I = 1; I <= AllocatedSize; I++)
    {
        Matrix->ExtToIntColMap[I] = -1;
        Matrix->ExtToIntRowMap[I] = -1;
    }
    Matrix->ExtToIntColMap[0] = 0;
    Matrix->ExtToIntRowMap[0] = 0;
#endif

    /* SRW */
#if BUILDHASH
    Matrix->ElementHashTab = NULL;
#endif
#if BITFIELD
    Matrix->BitField = NULL;
#endif

/* Allocate space for fill-ins and initial set of elements. */
    InitializeElementBlocks( Matrix, SPACE_FOR_ELEMENTS*AllocatedSize,
                                     SPACE_FOR_FILL_INS*AllocatedSize );
    if (Matrix->Error == spNO_MEMORY)
        goto MemoryError;

    return (char *)Matrix;

MemoryError:

/* Deallocate matrix and return no pointer to matrix if there is not enough
   memory. */
    *pError = spNO_MEMORY;
    spDestroy( (char *)Matrix);
    return NULL;
}









/*
 *  ELEMENT ALLOCATION
 *
 *  This routine allocates space for matrix elements. It requests large blocks
 *  of storage from the system and doles out individual elements as required.
 *  This technique, as opposed to allocating elements individually, tends to
 *  speed the allocation process.
 *
 *  >>> Returned:
 *  A pointer to an element.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *      A pointer to the first element in the group of elements being allocated.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

ElementPtr
spcGetElement( MatrixPtr Matrix )

{
ElementPtr  pElement;

/* Begin `spcGetElement'. */

/* Allocate block of MatrixElements if necessary. */
    if (Matrix->ElementsRemaining == 0)
    {
        pElement=0;
        spAlloc( Matrix, &pElement, sizeof(struct MatrixElement)*ELEMENTS_PER_ALLOCATION, AllocTypeMatrixElement  );
        if (Matrix->Error == spNO_MEMORY) return NULL;
        Matrix->ElementsRemaining = ELEMENTS_PER_ALLOCATION;
        Matrix->NextAvailElement = pElement;
    }

/* Update Element counter and return pointer to Element. */
    Matrix->ElementsRemaining--;
    return Matrix->NextAvailElement++;

}








/*
 *  ELEMENT ALLOCATION INITIALIZATION
 *
 *  This routine allocates space for matrix fill-ins and an initial set of
 *  elements.  Besides being faster than allocating space for elements one
 *  at a time, it tends to keep the fill-ins physically close to the other
 *  matrix elements in the computer memory.  This keeps virtual memory paging
 *  to a minimum.
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *  InitialNumberOfElements  <input> (int)
 *      This number is used as the size of the block of memory, in
 *      MatrixElements, reserved for elements. If more than this number of
 *      elements are generated, then more space is allocated later.
 *  NumberOfFillinsExpected  <input> (int)
 *      This number is used as the size of the block of memory, in
 *      MatrixElements, reserved for fill-ins. If more than this number of
 *      fill-ins are generated, then more space is allocated, but they may
 *      not be physically close in computer's memory.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *      A pointer to the first element in the group of elements being allocated.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

static void
InitializeElementBlocks( MatrixPtr Matrix, int InitialNumberOfElements,
                         int NumberOfFillinsExpected )
{
/* Begin `InitializeElementBlocks'. */

    /* Allocate block of MatrixElements for elements. */
    spAlloc( Matrix, &Matrix->NextAvailElement, sizeof(struct MatrixElement)*InitialNumberOfElements, AllocTypeMatrixElement );
    if (Matrix->Error == spNO_MEMORY) return;
    Matrix->ElementsRemaining = InitialNumberOfElements;

    /* Allocate block of MatrixElements for fill-ins. */
    spAlloc( Matrix, &Matrix->NextAvailFillin, sizeof(struct MatrixElement)*NumberOfFillinsExpected, AllocTypeMatrixElement );
    if (Matrix->Error == spNO_MEMORY) return;
    Matrix->FillinsRemaining = NumberOfFillinsExpected;

    /* Allocate a fill-in list structure. */
    spAlloc( Matrix, &Matrix->FirstFillinListNode, sizeof(struct FillinListNodeStruct), AllocTypeFillinListNodeStruct );
    if (Matrix->Error == spNO_MEMORY) return;
    Matrix->LastFillinListNode = Matrix->FirstFillinListNode;

    Matrix->FirstFillinListNode->pFillinList = Matrix->NextAvailFillin;
    Matrix->FirstFillinListNode->NumberOfFillinsInList =NumberOfFillinsExpected;
    Matrix->FirstFillinListNode->Next = NULL;

    return;
}










/*
 *  FILL-IN ALLOCATION
 *
 *  This routine allocates space for matrix fill-ins. It requests large blocks
 *  of storage from the system and doles out individual elements as required.
 *  This technique, as opposed to allocating elements individually, tends to
 *  speed the allocation process.
 *
 *  >>> Returned:
 *  A pointer to the fill-in.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

ElementPtr
spcGetFillin( MatrixPtr Matrix )

{
struct FillinListNodeStruct *pListNode;

/* Begin `spcGetFillin'. */
ElementPtr eptr;
#if NOT STRIP OR LINT
    if (Matrix->FillinsRemaining == 0)
        return spcGetElement( Matrix );
#endif
#if STRIP OR LINT

    if (Matrix->FillinsRemaining == 0)
    {   pListNode = Matrix->LastFillinListNode;

/* First see if there are any stripped fill-ins left. */
        if (pListNode->Next != NULL)
        {   Matrix->LastFillinListNode = pListNode = pListNode->Next;
            Matrix->FillinsRemaining = pListNode->NumberOfFillinsInList;
            Matrix->NextAvailFillin = pListNode->pFillinList;
        }
        else
        {
/* Allocate block of fill-ins. */
            spAlloc( Matrix, &Matrix->NextAvailFillin, sizeof(struct MatrixElement)*ELEMENTS_PER_ALLOCATION, AllocTypeMatrixElement );
            if (Matrix->Error == spNO_MEMORY) return NULL;
            Matrix->FillinsRemaining = ELEMENTS_PER_ALLOCATION;

/* Allocate a fill-in list structure. */
            spAlloc( Matrix, &pListNode->Next,sizeof(struct FillinListNodeStruct), AllocTypeFillinListNodeStruct  );
            if (Matrix->Error == spNO_MEMORY) return NULL;
            Matrix->LastFillinListNode = pListNode = pListNode->Next;

            pListNode->pFillinList = Matrix->NextAvailFillin;
            pListNode->NumberOfFillinsInList = ELEMENTS_PER_ALLOCATION;
            pListNode->Next = NULL;
        }
    }
#endif

/* Update Fill-in counter and return pointer to Fill-in. */
    ASSERT(Matrix->NextAvailFillin!=NULL)

    eptr = Matrix->NextAvailFillin;
    Matrix->FillinsRemaining--;
    if (Matrix->FillinsRemaining>0)
    {
      Matrix->NextAvailFillin++;
    }
    else
    {
      Matrix->NextAvailFillin=NULL;
    }

    return eptr;
}





/*
 *  RECORD A MEMORY ALLOCATION
 *
 *  This routine is used to record all memory allocations so that the memory
 *  can be freed later.
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *  AllocatedPtr  <input>  (char *)
 *      The pointer returned by malloc or calloc.  These pointers are saved in
 *      a list so that they can be easily freed.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */
void spAllocateBlockOfAllocationHashList( MatrixPtr Matrix )
{
    unsigned int  tcnt;
    spAllocationListHashPtr  ListPtr, ListPtr1;

    /* Begin `mulAllocateBlockOfAllocationList'. */
    ASSERT (Matrix!=NULL)

    /* Allocate block of records for allocation list. */
    ListPtr=calloc((spELEMENTS_PER_ALLOCATION+1), sizeof(struct spAllocationRecordHash));
    if (ListPtr == NULL)
    {

        return;
    }

    ListPtr1=ListPtr;
    /* Link all other new entries in EmptyAllocHashRecord chain */
    for (tcnt=0; tcnt < spELEMENTS_PER_ALLOCATION+1; tcnt++)
    {
         ListPtr->NextRecord = Matrix->EmptyAllocHashRecord;
         Matrix->EmptyAllocHashRecord= ListPtr;
         ListPtr++;
    }

    /* Register istself */
    spRecordAllocationFunc(Matrix, (char*) ListPtr1, (spELEMENTS_PER_ALLOCATION+1)*sizeof(struct spAllocationRecordHash), AllocTypeAllocationHashRecord, __FILE__, __LINE__ );
}

void spAllocateBlockOfAllocationList( MatrixPtr Matrix )
{
    unsigned int  tcnt;
    spAllocationListPtr  ListPtr, ListPtr1;

    /* Begin `mulAllocateBlockOfAllocationList'. */
    ASSERT (Matrix!=NULL)

    /* Allocate block of records for allocation list. */
    ListPtr=calloc((spELEMENTS_PER_ALLOCATION+1), sizeof(struct spAllocationRecord));
    if (ListPtr == NULL)
    {

        return;
    }

    ListPtr1=ListPtr;
    /* Link all other new entries in EmptyAllocHashRecord chain */
    for (tcnt=0; tcnt < spELEMENTS_PER_ALLOCATION+1; tcnt++)
    {
         ListPtr->NextRecord = Matrix->EmptyAllocRecord;
         Matrix->EmptyAllocRecord= ListPtr;
         ListPtr++;
    }

    /* Register itself */
    spRecordAllocationFunc(Matrix, (char*) ListPtr1, (spELEMENTS_PER_ALLOCATION+1)*sizeof(struct spAllocationRecord), AllocTypeAllocationRecord, __FILE__, __LINE__ );
}

void spRecordAllocationFunc( MatrixPtr Matrix, char *AllocatedPtr, unsigned long AllocatedSize, unsigned int AllocatedType, const char* filename, unsigned int linenumber )
{
    #define bittabsize MAX((1u<<spCalcHashWidth)/8,1)
    unsigned int  i, tsize, tsize1, done;
    unsigned char bittab[bittabsize];
    unsigned int thash;
    spAllocationListHashPtr tblockLast;
    spAllocationListPtr tallocblock;

    /*
     * If Allocated pointer is NULL, assume that malloc returned a NULL pointer,
    */
    ASSERT (Matrix!=NULL)

    ASSERT (AllocatedPtr!=NULL)
    ASSERT (AllocatedSize!=0)

    if ((AllocatedPtr == NULL) || (AllocatedSize==0))
    {
        return;
    }

    /* Replenish empty allocation record blocks if all are gone */
    if (Matrix->EmptyAllocRecord==NULL)
    {
       spAllocateBlockOfAllocationList( Matrix );
    }
    ASSERT (Matrix->EmptyAllocRecord!=NULL)

    /* Grab one block from Empty Alloc Record list */
    tallocblock=Matrix->EmptyAllocRecord;
    Matrix->EmptyAllocRecord=tallocblock->NextRecord;

    /* Chain in in Allocation List */
    if (Matrix->TopOfAllocationList!=NULL)
    {
       Matrix->TopOfAllocationList->PrevRecord=tallocblock;
    }
    tallocblock->NextRecord=Matrix->TopOfAllocationList;
    tallocblock->PrevRecord=NULL;
    Matrix->TopOfAllocationList=tallocblock;

    /* Fill with content */
    tallocblock->AllocatedPtr = AllocatedPtr;
    tallocblock->AllocatedSize = AllocatedSize;
    tallocblock->AllocatedType= AllocatedType;
    tallocblock->AllocatedDebugLineNumber=linenumber;
    tallocblock->AllocatedDebugFilename=filename;


    /* Create hash table entries and link them */
    tsize=AllocatedSize;
    for (i=0; i<bittabsize; i++)
    {
      bittab[i]=0;
    }
    tsize1=0; done=0;tblockLast=NULL;
    while (done==0)
    {
      unsigned int bitaddr, bitmask;

      /* Calculate hash for this address*/
      thash=spCalcHash( &AllocatedPtr[tsize1/sizeof(char)]  );

      /* Calculate bit position in bit map */
      bitaddr=thash>>3;
      bitmask=1<<(thash&7);

      /* If hash didn't occur so far add it */
      if (( bittab[ bitaddr] & bitmask) ==0)
      {
       spAllocationListHashPtr tblock;

       /* Replenish empty hash record blocks if all are gone */
       if (Matrix->EmptyAllocHashRecord==NULL)
       {
         spAllocateBlockOfAllocationHashList( Matrix );
       }
       ASSERT (Matrix->EmptyAllocHashRecord!=NULL)

       /* Remove from EmptyAllocHashRecord chain */
       tblock=Matrix->EmptyAllocHashRecord;
       Matrix->EmptyAllocHashRecord=tblock->NextRecord;

       /* Insert record in double linked list of specific hash chain */
       tblock->NextRecord = Matrix->TopOfAllocationHashList[thash];
       tblock->PrevRecord = NULL;
       if (Matrix->TopOfAllocationHashList[thash]!=NULL)
       {
          Matrix->TopOfAllocationHashList[thash]->PrevRecord=tblock;
       }
       Matrix->TopOfAllocationHashList[thash] = tblock;
       tblock->hashtable= thash;
       tblock->AllocRecord= tallocblock;

       /* Insert record in double linked list for identical entries in different hash tables */
       tblock->PrevPeerRecord=tblockLast;
       if (tblockLast!=NULL)
       {
         tblockLast->NextPeerRecord=tblock;
       }
       tblockLast=tblock;

       /* Set bit in bitmap to indicte this hash occured */
       bittab[ bitaddr ] |= bitmask ;
      }

      /* has changes every n addresses, so get to next address that changes hash */
      if (tsize1<tsize-1)
      {
        tsize1+=(1<<spCalcHashGridSize);
        if (tsize1>=tsize)
        {
          tsize1=tsize-1;
        }
      }
      else
      {
        done=1;
      }
    }
}

void spRemoveRecordAllocationFunc (MatrixPtr Matrix, char *AllocatedPtr, const char* filename, unsigned int linenumber)
{
   /* Begin `mulRemoveRecordAllocationFunc'. */
   ASSERT (Matrix!=NULL)
   spAllocationListHashPtr ListPtr;
   spAllocationListPtr tallocblock;

   if (AllocatedPtr == NULL)
    {
        return;
    }

    ListPtr=spFindRecordAllocation (Matrix, AllocatedPtr );
    ASSERT (ListPtr!=NULL)
    tallocblock=ListPtr->AllocRecord;

    /* Remove hash table entries */

    /* Move at the beginning of the peer hash entry chain */
    while (ListPtr->PrevPeerRecord!=NULL)
    {
      ListPtr=ListPtr->PrevPeerRecord;
    }

    /* Remove one by the other till nothing is left */
    while (ListPtr!=NULL)
    {
      spAllocationListHashPtr EListPtr;

      /* Just to be sure no memory curruption is ongoing */
      ASSERT(ListPtr->AllocRecord==tallocblock)

      /* Remove element from double link chained table */
      if (ListPtr->PrevRecord==NULL)
      {
        Matrix->TopOfAllocationHashList[ListPtr->hashtable]=ListPtr->NextRecord;
      }
      else
      {
        ListPtr->PrevRecord->NextRecord=ListPtr->NextRecord;
      }
      if (ListPtr->NextRecord!=NULL)
      {
        ListPtr->NextRecord->PrevRecord=ListPtr->PrevRecord;
      }

      /* Progress with next element */
      EListPtr=ListPtr;
      ListPtr=ListPtr->NextPeerRecord;

      /* Empty element */
      memset(EListPtr, 0, sizeof(struct spAllocationRecordHash));

      /* Recycle entry */
      EListPtr->NextRecord=Matrix->EmptyAllocHashRecord;
      Matrix->EmptyAllocHashRecord=EListPtr;
    }

    /* Remove allocation entry from chain*/
    if (tallocblock->PrevRecord==NULL)
    {
      Matrix->TopOfAllocationList=tallocblock->NextRecord;
    }
    else
    {
      tallocblock->PrevRecord->NextRecord= tallocblock->NextRecord;
    }
    if (tallocblock->NextRecord!=NULL)
    {
      tallocblock->NextRecord->PrevRecord=tallocblock->PrevRecord;
    }

    /* Empty block */
    memset(tallocblock, 0, sizeof(struct spAllocationRecord));

    /* Link it to empty chain */
    tallocblock->NextRecord=Matrix->EmptyAllocRecord;
    Matrix->EmptyAllocRecord=tallocblock;
}




spAllocationListHashPtr spFindRecordAllocation (MatrixPtr Matrix, char *Ptr )
{
/* Begin `mulFindRecordAllocation'. */

  spAllocationListHashPtr  ListPtr;
  ASSERT (Matrix!=NULL)
  unsigned int thash;
  unsigned int cnt;

    if (Ptr==NULL)
    {
        return NULL;
    }

    thash=spCalcHash(Ptr);
    ListPtr = Matrix->TopOfAllocationHashList[thash];
    while (ListPtr != NULL)
    {
        ASSERT(ListPtr->AllocRecord!=NULL)
        ASSERT(ListPtr->AllocRecord->AllocatedPtr!=NULL)
        ASSERT(ListPtr->AllocRecord->AllocatedSize>0)

        if (
            ( ListPtr->AllocRecord->AllocatedPtr                                                    <=Ptr) &&
            (&ListPtr->AllocRecord->AllocatedPtr[ListPtr->AllocRecord->AllocatedSize/sizeof(char)]  > Ptr)
           )
        {
            return ListPtr;
        }
        else
        {
          ListPtr = ListPtr->NextRecord;
        }
    }

    return NULL;
}

void spAllocFunc( MatrixPtr Matrix, void** AllocatedPtr, unsigned long AllocatedSize, unsigned int AllocatedType, const char* filename, unsigned int linenumber )
{

  ASSERT (Matrix!=NULL)

  if ((*AllocatedPtr)!=NULL)
  {
    free(*AllocatedPtr);
    spRemoveRecordAllocationFunc(Matrix, *AllocatedPtr, filename, linenumber);
  }

  *AllocatedPtr=NULL;

  if (AllocatedSize==0)
  {
    return;
  }

  *AllocatedPtr=calloc(AllocatedSize/sizeof(char),sizeof(char));
  spRecordAllocationFunc(Matrix,(char*)*AllocatedPtr,AllocatedSize,AllocatedType,filename,linenumber);
}

void spReallocFunc (MatrixPtr Matrix, void** pPtr, unsigned long Size, const char* filename, unsigned int linenumber )
{
    void *oldPtr;
    spAllocationListHashPtr  ListPtr;
    unsigned int typ;

   ASSERT (Matrix!=NULL)
   ASSERT (pPtr!=NULL)
   ASSERT (*pPtr!=NULL)
   ASSERT (Size>0)

   oldPtr= *pPtr;

    ListPtr = spFindRecordAllocation(Matrix, oldPtr);
    ASSERT (ListPtr!=NULL)
    typ= ListPtr->AllocRecord->AllocatedType;

    (*pPtr) = realloc(*pPtr, Size);

    spRemoveRecordAllocationFunc(Matrix, (char*)oldPtr, filename, linenumber);

    spRecordAllocationFunc(Matrix, (char*)*pPtr, Size, typ, filename, linenumber);

}

void spFreeFunc (MatrixPtr Matrix, void** pPtr, const char* filename, unsigned int linenumber )
{
   spAllocationListHashPtr ListPtr;
   spAllocationListHashPtr tblock;
   unsigned int  i, tsize, tsize1, done;
   unsigned char bittab[(1u<<spCalcHashWidth)/8];
   char* ptr;

   ASSERT (Matrix!=NULL)
   ASSERT (pPtr!=NULL)
   ASSERT (*pPtr!=NULL)

   ListPtr=spFindRecordAllocation(Matrix, *pPtr);
   ASSERT (ListPtr!=NULL)

   ptr=(char*)*pPtr;

   /* Free memory block */
   free(*pPtr);

   /* Delete link to block */
   (*pPtr) = NULL;

    spRemoveRecordAllocationFunc (Matrix, ptr,filename,linenumber);
}

unsigned long spGetMemoryBlockLength (MatrixPtr Matrix, char* AllocatedPtr)
{
   spAllocationListHashPtr  ListPtr;
    ASSERT (Matrix!=NULL)

    if (AllocatedPtr == NULL)
    {
        return 0;
    }

    ListPtr = spFindRecordAllocation(Matrix, AllocatedPtr);
    if (ListPtr!=NULL)
    {
        return ListPtr->AllocRecord->AllocatedSize;
    }

    return 0;
}

void spDestroy( char *eMatrix )
{
  MatrixPtr Matrix = (MatrixPtr)eMatrix;

  spAllocationListPtr  ListPtr;
  unsigned int i;


  if (Matrix==NULL)
  {
    /* Nothing to do */
    return;
  }
  ASSERT( IS_SPARSE( Matrix ) );

    /* SRW */
#if BUILDHASH
    sph_destroy(Matrix);
#endif
#if BITFIELD
    ba_destroy(Matrix);
#endif

  /* Make a copy of TopAllocationList entry */
  /* As soon block containing Matrix is gone */
  /* there is no access anymore */
  ListPtr=Matrix->TopOfAllocationList;

  /* Sequentially step through the list of allocated pointers freeing pointers
   * along the way. */

  /* Run twice over list first run remove all data blocks */
  /* second run remove allocation records itself */
  for (i=0; i<2; i++)
  {
      unsigned int cnt;
      spAllocationListPtr  tListPtr;

      cnt=0;
      tListPtr = ListPtr;
      while (tListPtr != NULL)
      {
          if ((i==1) || (tListPtr->AllocatedType!=AllocTypeAllocationRecord))
          {
            void* tptr, *tptr1;

            /* Save pointer to block to release */
            tptr=tListPtr->AllocatedPtr;

            /* Remove element from chain */
            if (tListPtr->PrevRecord==NULL)
            {
               ListPtr=tListPtr->NextRecord;
            }
            else
            {
              tListPtr->PrevRecord->NextRecord=tListPtr->NextRecord;
            }

            if (tListPtr->NextRecord!=NULL)
            {
               tListPtr->NextRecord->PrevRecord=tListPtr->PrevRecord;
            }
            tptr1=(void*)tListPtr;
            tListPtr=tListPtr->NextRecord;

            /* For sanity clear entry */
            memset(tptr1, 0, sizeof(struct spAllocationRecord));

            /* If this block is a allocation record block */
            /* free releases it and you aren't allowed to access */
            /* beyond this point */
            free(tptr);
            cnt++;
          }
          else
          {
            tListPtr=tListPtr->NextRecord;
          }
     }
//     printf ("%u memory blocks released \n",cnt);
  }
  return;
}

unsigned int spCalcHash(void* ptr)
{

  unsigned long addr;
  unsigned int thash;

  addr=((unsigned long)ptr)>>spCalcHashGridSize;

  thash = ((addr                      )&((1u<<spCalcHashWidth)-1));
  thash^= ((addr>>(  spCalcHashWidth))&((1u<<spCalcHashWidth)-1));
  thash^= ((addr>>(2*spCalcHashWidth))&((1u<<spCalcHashWidth)-1));
  thash^= ((addr>>(3*spCalcHashWidth))&((1u<<spCalcHashWidth)-1));
  thash^= ((addr>>(4*spCalcHashWidth))&((1u<<spCalcHashWidth)-1));
  thash^= ((addr>>(5*spCalcHashWidth))&((1u<<spCalcHashWidth)-1));
  thash^= ((addr>>(6*spCalcHashWidth))&((1u<<spCalcHashWidth)-1));

  return thash;
}





/*
 *  RETURN MATRIX ERROR STATUS
 *
 *  This function is used to determine the error status of the given matrix.
 *
 *  >>> Returned:
 *      The error status of the given matrix.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      The matrix for which the error status is desired.
 */

int
spError( char *eMatrix )

{
/* Begin `spError'. */

    if (eMatrix != NULL)
    {   ASSERT(((MatrixPtr)eMatrix)->ID == SPARSE_ID);
        return ((MatrixPtr)eMatrix)->Error;
    }
    else return spNO_MEMORY;   /* This error may actually be spPANIC,
                                * no way to tell. */
}









/*
 *  WHERE IS MATRIX SINGULAR
 *
 *  This function returns the row and column number where the matrix was
 *  detected as singular or where a zero was detected on the diagonal.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      The matrix for which the error status is desired.
 *  pRow  <output>  (int *)
 *      The row number.
 *  pCol  <output>  (int *)
 *      The column number.
 */

void
spWhereSingular( char *eMatrix, int *pRow, int *pCol )

{
MatrixPtr Matrix = (MatrixPtr)eMatrix;

/* Begin `spWhereSingular'. */
    ASSERT( IS_SPARSE( Matrix ) );

    if (Matrix->Error == spSINGULAR OR Matrix->Error == spZERO_DIAG)
    {   *pRow = Matrix->SingularRow;
        *pCol = Matrix->SingularCol;
    }
    else *pRow = *pCol = 0;
    return;
}






/*
 *  MATRIX SIZE
 *
 *  Returns the size of the matrix.  Either the internal or external size of
 *  the matrix is returned.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to matrix.
 *  External  <input>  (BOOLEAN)
 *      If External is set true, the external size , i.e., the value of the
 *      largest external row or column number encountered is returned.
 *      Otherwise the true size of the matrix is returned.  These two sizes
 *      may differ if the TRANSLATE option is set true.
 */

int
spGetSize( char *eMatrix, BOOLEAN External )

{
MatrixPtr Matrix = (MatrixPtr)eMatrix;

/* Begin `spGetSize'. */
    ASSERT( IS_SPARSE( Matrix ) );

#if TRANSLATE
    if (External)
        return Matrix->ExtSize;
    else
        return Matrix->Size;
#else
    return Matrix->Size;
#endif
}








/*
 *  SET MATRIX COMPLEX OR REAL
 *
 *  Forces matrix to be either real or complex.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to matrix.
 */

void
spSetReal( char *eMatrix )

{
/* Begin `spSetReal'. */

    ASSERT( IS_SPARSE( (MatrixPtr)eMatrix ) AND REAL);
    ((MatrixPtr)eMatrix)->Complex = NO;
    return;
}


void
spSetComplex( char *eMatrix )

{
/* Begin `spSetComplex'. */

    ASSERT( IS_SPARSE( (MatrixPtr)eMatrix ) AND spCOMPLEX);
    ((MatrixPtr)eMatrix)->Complex = YES;
    return;
}









/*
 *  ELEMENT OR FILL-IN COUNT
 *
 *  Two functions used to return simple statistics.  Either the number
 *  of total elements, or the number of fill-ins can be returned.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to matrix.
 */

int
spFillinCount( char *eMatrix )

{
/* Begin `spFillinCount'. */

    ASSERT( IS_SPARSE( (MatrixPtr)eMatrix ) );
    return ((MatrixPtr)eMatrix)->Fillins;
}


int
spElementCount( char *eMatrix )

{
/* Begin `spElementCount'. */

    ASSERT( IS_SPARSE( (MatrixPtr)eMatrix ) );
    return ((MatrixPtr)eMatrix)->Elements;
}

typedef struct memtable_t {
	void* old;
	void* new;
	unsigned long length;
        unsigned int type;
	} memtable_t;

#define UpdatePointer(ptr,memtable,memtableentries,hashtranstable) UpdatePointerFunc((void**)ptr,memtable, memtableentries, hashtranstable, __FILE__,__LINE__)
unsigned int UpdatePointerFunc(void** ptr, memtable_t* memtable, unsigned int memtableentries, unsigned long* hashtranstable, const char* filename, unsigned int linenumber)
{
  unsigned long lptr, offs;
  unsigned int i,j;
  unsigned int thash;
  unsigned int elems;

  if (*ptr==NULL)
  {
    return 0;
  }

  lptr=(unsigned long)*ptr;

  thash=spCalcHash(*ptr);
  offs=hashtranstable[thash];
  elems=hashtranstable[thash+1]-offs;

  for (j=0; j<elems; j++)
  {
    /* Hash transformation table lists per hash value the linked allocation entries */
    /* If hash works well (equal distribution) speed up goes with 2+^(hash width), e.g. 9 --> 512 hash values, so */
    /* there is best case only 1/512 of the overall allocation list to be scanned. MASSIVE SPEED UP ** YUMMY ** */
    i=hashtranstable[offs+j];

    if ((lptr>=(unsigned long)memtable[i].old) && (lptr<(unsigned long)memtable[i].old+memtable[i].length))
    {
       lptr=lptr-(unsigned long)memtable[i].old+(unsigned long)memtable[i].new;
       *ptr=(void*)lptr;
       return 0;
    }

  }
  fprintf(stderr,"Pointer translation not successful %s:%u!\n",filename,linenumber);
  return 1;

/* run over all elements to find the maching one  - that takes quite a while especially with big data

  unsigned long lptr;
  unsigned int i;

  if (*ptr==NULL)
  {
    return 0;
  }
  lptr=(unsigned long)*ptr;

  for (i=0;i<memtableentries; i++)
  {
    if ((lptr>=(unsigned long)memtable->old) && (lptr<(unsigned long)memtable->old+memtable->length))
    {
       lptr=lptr-(unsigned long)memtable->old+(unsigned long)memtable->new;
       *ptr=(void*)lptr;
       return 0;
    }
    memtable++;
  }

  #ifdef MPI
  {
    int MPIrank;
    MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);
    fprintf(stderr, "Rank %u: Pointer translation 0x%08lx not successful %s:%u!\n",MPIrank, lptr, filename,linenumber);
   }
  #else
   fprintf(stderr, "Pointer translation 0x%08lx not successful %s:%u!\n", lptr, filename,linenumber);
  #endif
  return 1;*/
}

unsigned int spMPIBlockNmbr=1;

/* */
/* Function to redirect data stream for load/save from/to FILE or MPI */
/* */
void spComm(unsigned int mode, FILE* fop, void* ptr, unsigned long size)
{


  switch (mode)
  {
    #ifdef MPI
    case spCommBcastMPI: // MPI Bcast
     MPI_Bcast(ptr, size/sizeof(char), MPI_CHAR, 0, *((MPI_Comm*)fop));
     break;

    case spCommtoMPI: // MPI Send
//     printf("MPI_Send spMPIBlockNmbr %u size %u to %u \n", spMPIBlockNmbr, size, *((int*)fop)); fflush (stdout);
     MPI_Send(ptr, size/sizeof(char), MPI_CHAR,
               *((int*)fop),  (spMPIBlockNmbr++), MPI_COMM_WORLD);
     break;

    case spCommfromMPI: // MPI Receive
    {
      MPI_Status status;
//    printf("MPI_Recv spMPIBlockNmbr %u size %u from %u \n", spMPIBlockNmbr, size, *((int*)fop)); fflush (stdout);
     MPI_Recv(ptr, size/sizeof(char), MPI_CHAR,
               *((int*)fop),  (spMPIBlockNmbr++), MPI_COMM_WORLD, &status);
    }
     break;

    #endif

   case spCommtoFile: // File save
     fwrite(ptr, size, 1, fop);
     break;

   case spCommfromFile: // File read
     fread(ptr, size, 1, fop);
     break;

   default:
     fprintf(stderr, "spMatrixComm: Unknown communication target");
     exit (-1);
     break;
  }
 /**/
  if (spMPIBlockNmbr>=4094)
  {
    spMPIBlockNmbr=0;
  }
}

/* Load / Save / MPI transfer function */
/* mode: spCommtoMPI Copy over MPI */
/*       spCommtoFile Save */
/*       spCommfromFile Load */

void spMatrixComm (char** ppsparMatrix, unsigned int mode, FILE* fop)
{
   int MPIrank;
   unsigned long i;
   MatrixPtr  Matrix;
   spAllocationListPtr  ListPtr;
   memtable_t *memtable;
   unsigned int memtableentries;

   ASSERT(ppsparMatrix!=NULL)

   switch (mode)
   {
   case spCommBcastMPI: // MPI Bcast
     #ifdef MPI
       ASSERT (fop!=NULL)
       MPI_Comm_rank (*((MPI_Comm*)fop), &MPIrank);	/* get current process id */
     #else
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     break;

   case spCommtoMPI:
     #ifndef MPI
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     spMPIBlockNmbr=1;
   case spCommtoFile: // File save
     MPIrank=0;
     break;

   case spCommfromMPI:
     #ifndef MPI
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     spMPIBlockNmbr=1;
   case spCommfromFile: // File read
     MPIrank=1;
     break;

   default:
     fprintf(stderr, "spMatrixComm: Unknown communication target");
     exit (-1);
     break;
   }

   if (MPIrank>0)
   {
     spDestroy(*ppsparMatrix);
     *ppsparMatrix=NULL;
   }

//  printf ("%u: Entering Matrix 0x%08lx\n", MPIrank, (unsigned long)(*ppsparMatrix) ); fflush(stdout);
   /* Copy content */
   spComm(mode, fop,
          (void*)(ppsparMatrix),
          sizeof(MatrixPtr*));

   if ((*ppsparMatrix)==NULL)
   {
     return;
   }

//   printf ("%u: Counting elements Matrix 0x%08lx\n", MPIrank, (unsigned long)(*ppsparMatrix) ); fflush(stdout);
   {
      MatrixPtr Matrix;
      spAllocationListPtr  ListPtr;
      spAllocationListHashPtr HashListPtr;
      memtable_t *memtable;
      unsigned long memtableentries, memtable_size;
      unsigned long hashtableentries[1<<spCalcHashWidth],hashtableentries_all;
      unsigned long *hashtranstable,hashtranstable_size;

      if (MPIrank==0)
      {
        unsigned int cnt;

        Matrix = *((MatrixPtr*)ppsparMatrix);
        ListPtr = Matrix->TopOfAllocationList;
        memtableentries=0;
        while (ListPtr != NULL)
        {
            /* Number them through for Hash table transfer */
            ListPtr->AllocNumber=memtableentries;
            /* Work on next element */
            ListPtr = ListPtr->NextRecord;
            memtableentries++;
        }

        hashtableentries_all=0;
        for (cnt=0; cnt< (1<<spCalcHashWidth); cnt++ )
        {
           unsigned long hentries;
           spAllocationListHashPtr HashListPtr;

           hentries=0;
           HashListPtr=Matrix->TopOfAllocationHashList[cnt];
           while (HashListPtr!=NULL)
           {
             hentries++;
             HashListPtr=HashListPtr->NextRecord;
           }
           hashtableentries[cnt]=hentries;
           hashtableentries_all+=hentries;
        }

      }
      else
      {
        /* Generate a new main data element = allocation function links hash and allocation records in there */
        /* Do not register it: As soon the received main data structure gets processed those links get copied and this structure released again */
        *ppsparMatrix = (char*)calloc(1, sizeof(struct MatrixFrame));
      }
      Matrix = *((MatrixPtr*)ppsparMatrix);

      spComm(mode, fop, (void*)(&memtableentries),    sizeof(unsigned long));
      spComm(mode, fop, (void*)&hashtableentries_all, sizeof(unsigned long));


      memtable_size=(unsigned long)sizeof(struct memtable_t)*((unsigned long)memtableentries);
      memtable=malloc(memtable_size);
      ASSERT (memtable!=NULL);

      if (MPIrank==0)
      {
	    memtable_t *tmemtable;

        ListPtr = Matrix->TopOfAllocationList;
	    tmemtable=memtable;
        while (ListPtr!=NULL)
        {
          tmemtable->old=ListPtr->AllocatedPtr;
          tmemtable->new=0;
   	      tmemtable->length=ListPtr->AllocatedSize;
          tmemtable->type=ListPtr->AllocatedType;
          ListPtr = ListPtr->NextRecord;
          tmemtable++;
        }
      }

      hashtranstable_size=(unsigned long)sizeof(unsigned long) * (hashtableentries_all+(1u<<spCalcHashWidth)+1);
      hashtranstable=malloc(hashtranstable_size);
      ASSERT (hashtranstable!=NULL)
      if (MPIrank==0)
      {
         unsigned int cnt;
         unsigned long entry;
         unsigned long *tlongptr;
         unsigned long tlcnt;
         entry=(1u<<spCalcHashWidth)+1;
         tlongptr=hashtranstable;
         tlcnt=0;
         for (cnt=0; cnt< (1u<<spCalcHashWidth); cnt++)
         {
            *(tlongptr++)=entry;
            entry+=hashtableentries[cnt];
            tlcnt++;
         }
         *(tlongptr++)=entry;
         tlcnt++;

         for (cnt=0; cnt< (1u<<spCalcHashWidth); cnt++)
         {
           spAllocationListHashPtr HashListPtr;

           // Assure consistency of table
           ASSERT(tlcnt==hashtranstable[cnt]);
           HashListPtr=Matrix->TopOfAllocationHashList[cnt];
           while (HashListPtr!=NULL)
           {
             *(tlongptr++)=HashListPtr->AllocRecord->AllocNumber;
             tlcnt++;
             HashListPtr=HashListPtr->NextRecord;
           }
         }
         // Assure consistency of table
         ASSERT(tlcnt==hashtranstable[cnt])
      }

      spComm(mode, fop, (void*)memtable, memtable_size);
//      printf ("%u: MemTable copied\n", MPIrank); fflush(stdout);
      spComm(mode, fop, (void*)hashtranstable, hashtranstable_size);

      if (MPIrank>0)
      {
	    unsigned int memsize;
	    memtable_t *tmemtable;
        memsize=0;
        tmemtable=memtable;
        for (i=0; i<memtableentries; i++)
        {
           if ((tmemtable->type!=AllocTypeAllocationHashRecord) &&
               (tmemtable->type!=AllocTypeAllocationRecord))
           {
               ASSERT (tmemtable->length>0)

               /* Allocate memory */
               tmemtable->new=malloc(tmemtable->length);
               /* Register record and built immediatly local hash structure */
               spRecordAllocationFunc( Matrix,
                                       (char*)tmemtable->new,
                                        tmemtable->length,
                                        tmemtable->type,
                                        __FILE__,__LINE__);
               memsize+=tmemtable->length;
           }
           tmemtable++;
        }
//        printf ("%u: Memory allocated (%u kB)\n", MPIrank, memsize>>10); fflush(stdout);
      }
      {
        #ifdef STATISTIC
        unsigned long stat_generic;
        unsigned long stat_matrixframe;
        unsigned long stat_matrixelement;
        unsigned long stat_fillinlistnodestruc;
        unsigned long stat_allocrcord;
        unsigned long stat_ptrlist;
        #endif
    	memtable_t *tmemtable;

        #ifdef STATISTIC
        stat_generic=0;stat_matrixframe=0;stat_matrixelement=0;stat_fillinlistnodestruc=0;stat_allocrcord=0;stat_ptrlist=0;
        #endif

        tmemtable=memtable;
        for (i=0; i<memtableentries; i++)
        {
           #ifdef MPI
           if (mode==spCommBcastMPI)
           {
             // Dirty: openmpi-1.7.3-1.fc20
             // For what ever reason more than 8192 Bcasts aren't allowed to be pending ...
             if ((i&2047)==0)
             {
               MPI_Barrier (*(MPI_Comm*)fop );
             }
           }
          #endif
          if ((tmemtable->type!=AllocTypeAllocationHashRecord) &&
               (tmemtable->type!=AllocTypeAllocationRecord))
           {
             void* temp;
             unsigned int err;

             ASSERT ((tmemtable->length>0))
             if (MPIrank>0)
             {
                temp=tmemtable->new;
             }
             else
             {
               temp=tmemtable->old;
             }
             spComm(mode, fop, (void*)temp, tmemtable->length);

             if (MPIrank>0)
             {
             err=0;
	         /* Update pointers in copied content */
             switch (tmemtable->type)
             {
             case AllocTypeGeneric:
               #ifdef STATISTIC
                 stat_generic++;
               #endif
               /* Nothing to do */
               break;

             case AllocTypePtrArray:
               {
                 void** tptr;
                 unsigned long tcount;
                 unsigned int i;
                 #ifdef STATISTIC
                   stat_ptrlist++;
                 #endif
                 tptr=(void**)temp;
                 tcount=tmemtable->length/sizeof(void*);
                 for (i=0; i<tcount; i++)
                 {
                    UpdatePointer(&(*tptr),memtable,memtableentries,hashtranstable);
                    tptr++;
                 }
               }
               break;

             case AllocTypeMatrixFrame:
               {
                 MatrixPtr tMatrix;
                 unsigned int i;

                 #ifdef STATISTIC
                   stat_matrixframe++;
                 #endif
                 tMatrix=(MatrixPtr)temp;
                 ASSERT(tmemtable->length==sizeof(struct MatrixFrame))
                 UpdatePointer(&tMatrix->Intermediate,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->Diag,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->FirstInCol,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->FirstInRow,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->DoCmplxDirect,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->DoRealDirect,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->ExtToIntColMap,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->ExtToIntRowMap,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->IntToExtColMap,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->IntToExtRowMap,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->MarkowitzRow,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->MarkowitzCol,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->MarkowitzProd,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->NextAvailElement,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->NextAvailFillin,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->FirstFillinListNode,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->LastFillinListNode,memtable,memtableentries,hashtranstable);
#if BUILDHASH
                 UpdatePointer(&tMatrix->ElementHashTab,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&tMatrix->he_blocks,memtable,memtableentries,hashtranstable);
#endif // BUILDHASH
#if BITFIELD
                 UpdatePointer(&tMatrix->BitField,memtable,memtableentries,hashtranstable);
#endif

                 /* Copy memory management links */
                 for (i=0; i<(1u<<spCalcHashWidth); i++)
                 {
                   tMatrix->TopOfAllocationHashList[i]=Matrix->TopOfAllocationHashList[i];
                 }
                 tMatrix->TopOfAllocationList=Matrix->TopOfAllocationList;
                 tMatrix->EmptyAllocHashRecord=Matrix->EmptyAllocHashRecord;
                 tMatrix->EmptyAllocRecord=Matrix->EmptyAllocRecord;

                 free (Matrix);

                 (*ppsparMatrix)=(char*)tMatrix;
               }
               break;

#if BUILDHASH

            case AllocTypespHtab:
            {
                struct spHtab *tab;
                unsigned int i, tcount;

                 tcount=tmemtable->length/sizeof(struct spHtab);
                 tab=(struct spHtab *)temp;

                 for (i=0; i<tcount; i++)
                 {
                   UpdatePointer(&tab->entries,memtable,memtableentries,hashtranstable);
                   tab++;
                 }
            }
            break;

            case AllocTypespHelt:
            {
                struct spHelt *elt;
                unsigned int i, tcount;

                 tcount=tmemtable->length/sizeof(struct spHelt);
                 elt=(struct spHelt *)temp;

                 for (i=0; i<tcount; i++)
                 {
                   UpdatePointer(&elt->eptr,memtable,memtableentries,hashtranstable);
                   UpdatePointer(&elt->next,memtable,memtableentries,hashtranstable);
                   elt++;
                 }
            }
            break;

            case AllocTypespHeblk:
            {
                struct heblk *eblk;
                unsigned int i, tcount;

                 tcount=tmemtable->length/sizeof(struct heblk);
                 eblk=(struct heblk *)temp;

                 for (i=0; i<tcount; i++)
                 {
                   unsigned int j;

                   for (j=0; j<HEBLKSZ; j++)
                   {
                      UpdatePointer(&eblk->elts[j].next,memtable,memtableentries,hashtranstable);
                      UpdatePointer(&eblk->elts[j].eptr,memtable,memtableentries,hashtranstable);
                    }
                   UpdatePointer(&eblk->next,memtable,memtableentries,hashtranstable);
                   eblk++;
                 }
            }
            break;


#endif // BUILDHASH

             case AllocTypeMatrixElement:
               {
                 ElementPtr  pElement;
                 unsigned int i;
                 unsigned int tcount;
                 pElement=(ElementPtr)temp;
                 tcount=tmemtable->length/sizeof(struct MatrixElement);
                 #ifdef STATISTIC
                   stat_matrixelement++;
                 #endif

                 for (i=0; i<tcount; i++)
                 {
                   UpdatePointer(&pElement->NextInRow,memtable,memtableentries,hashtranstable);
                   UpdatePointer(&pElement->NextInCol,memtable,memtableentries,hashtranstable);
                   #if INITIALIZE
                     UpdatePointer(&pElement->pInitInfo,memtable,memtableentries,hashtranstable);
                   #endif
                   pElement++;
                 }
               }
               break;

             case AllocTypeFillinListNodeStruct:
               {
                 struct FillinListNodeStruct  *pElement;
                 #ifdef STATISTIC
                   stat_fillinlistnodestruc++;
                 #endif
                 pElement=(struct FillinListNodeStruct  *)temp;
                 ASSERT(tmemtable->length==sizeof(struct FillinListNodeStruct))
                 UpdatePointer(&pElement->pFillinList,memtable,memtableentries,hashtranstable);
                 UpdatePointer(&pElement->Next,memtable,memtableentries,hashtranstable);
               }
               break;

               default:
                 fprintf(stderr,"CopysparMatrixoverMPI: Unknown AllocatedType %u\n",tmemtable->type);
                 exit(-1);
               }
             }
          }
          tmemtable++;
        }
#ifdef STATISTIC
        printf ("%u: MemTable linked content copied \n", MPIrank); fflush(stdout);
        printf ("%u: Generic Blocks %u\n", MPIrank, stat_generic);
        printf ("%u: Pointer lists %u\n", MPIrank, stat_ptrlist);
        printf ("%u: Matrix frame %u\n", MPIrank, stat_matrixframe);
        printf ("%u: Matrix element %u\n", MPIrank, stat_matrixelement);
        printf ("%u: fillinlistnodestruc %u\n", MPIrank,stat_fillinlistnodestruc);
        printf ("%u: AllocationRecord %u\n", MPIrank, stat_allocrcord);
#endif
     }
     free (memtable);
     free(hashtranstable);
   }
}

#define MarkAlloc(Matrix, ptr, flag) MarkAllocFunc(Matrix, (void*)ptr, flag, __FILE__, __LINE__)
int MarkAllocFunc (MatrixPtr Matrix, void* ptr, unsigned int flag, const char* filename, unsigned int linenumber)
{
    ASSERT (Matrix!=NULL)

    spAllocationListHashPtr  ListPtr;

    if (ptr==NULL)
      return 0;

    ListPtr = spFindRecordAllocation(Matrix, (char*)ptr);
    if (ListPtr!=NULL)
    {
      if (flag==0)
      {
        return 0;
      }

      /* Mark this entry */
      if (ListPtr->AllocRecord->Mark==0)
      {
        ListPtr->AllocRecord->Mark=1;
      }

      return 0;
    }
    fprintf(stderr, "Pointer (0x%08lx) not found in allocation list %s:%u\n", (unsigned long)ptr, filename, linenumber);
    return 1;
}

int spAnalyzeAllocatedMemory (char* psparMatrix)
{
    MatrixPtr Matrix;
    spAllocationListPtr  ListPtr;
    unsigned int found;
    unsigned int allocs;
    unsigned long allocsize;
    unsigned int MPIrank;

    if (psparMatrix==NULL)
    {
      return 0;
    }

    MPIrank=0;
    #ifdef MPI
    MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);
    #endif

//   printf ("%u: Counting elements Matrix 0x%08lx\n", MPIrank, (unsigned long)(psparMatrix) ); fflush(stdout);

    /* Clear all marks */
    allocs=0;
    Matrix=(MatrixPtr)psparMatrix;
    ListPtr=Matrix->TopOfAllocationList;
    while (ListPtr!=NULL)
    {
        ListPtr->Mark=0;
        allocs++;
        ListPtr=ListPtr->NextRecord;
    }

   printf ("%u: %u elements\n", MPIrank, allocs); fflush(stdout);

    /* Mark root structure and develop from there */
    MarkAlloc(Matrix, Matrix, 1);

    /* Run iterative over list until there isn't anything to mark anymore */
    found=1;
    while (found>0)
    {
      found=0;
      ListPtr=Matrix->TopOfAllocationList;
      while (ListPtr!=NULL)
      {
        if (ListPtr->Mark==1)
        {
            unsigned int err;
            found=1;
            ListPtr->Mark=2;
            err=0;
            switch (ListPtr->AllocatedType)
            {
            case AllocTypeGeneric:
              /* Nothing to do */
              break;

            case AllocTypePtrArray:
              {
                void** ptr;
                unsigned long tcount;
                unsigned int i;
                ptr=(void**)ListPtr->AllocatedPtr;
                tcount=ListPtr->AllocatedSize/sizeof(void*);
                for (i=0; i<tcount; i++)
                {
                   err|=MarkAlloc(Matrix,(*ptr),1);
                   ptr++;
                }
              }
              break;

               case AllocTypeMatrixFrame:
                 {
                   unsigned int i;
                   Matrix=(MatrixPtr)ListPtr->AllocatedPtr;
                   ASSERT(ListPtr->AllocatedSize==sizeof(struct MatrixFrame))
                   err|=MarkAlloc(Matrix,Matrix->Intermediate,1);
                   err|=MarkAlloc(Matrix,Matrix->Diag,1);
                   err|=MarkAlloc(Matrix,Matrix->FirstInCol,1);
                   err|=MarkAlloc(Matrix,Matrix->FirstInRow,1);
                   err|=MarkAlloc(Matrix,Matrix->DoCmplxDirect,1);
                   err|=MarkAlloc(Matrix,Matrix->DoRealDirect,1);
                   err|=MarkAlloc(Matrix,Matrix->ExtToIntColMap,1);
                   err|=MarkAlloc(Matrix,Matrix->ExtToIntRowMap,1);
                   err|=MarkAlloc(Matrix,Matrix->IntToExtColMap,1);
                   err|=MarkAlloc(Matrix,Matrix->IntToExtRowMap,1);
                   err|=MarkAlloc(Matrix,Matrix->MarkowitzRow,1);
                   err|=MarkAlloc(Matrix,Matrix->MarkowitzCol,1);
                   err|=MarkAlloc(Matrix,Matrix->MarkowitzProd,1);
                   err|=MarkAlloc(Matrix,Matrix->NextAvailElement,1);
                   err|=MarkAlloc(Matrix,Matrix->NextAvailFillin,1);
                   err|=MarkAlloc(Matrix,Matrix->FirstFillinListNode,1);
                   err|=MarkAlloc(Matrix,Matrix->LastFillinListNode,1);
                   for (i=0; i<(1u<<spCalcHashWidth); i++)
                   {
                     err|=MarkAlloc(Matrix,Matrix->TopOfAllocationHashList[i],1);
                   }
                   err|=MarkAlloc(Matrix,Matrix->EmptyAllocHashRecord,1);
                   err|=MarkAlloc(Matrix,Matrix->TopOfAllocationList,1);
                   err|=MarkAlloc(Matrix,Matrix->EmptyAllocRecord,1);
#if BUILDHASH
                 err|=MarkAlloc(Matrix,Matrix->ElementHashTab,1);
                 err|=MarkAlloc(Matrix,Matrix->he_blocks,1);
#endif // BUILDHASH
#if BITFIELD
                 err|=MarkAlloc(Matrix,Matrix->BitField,1);
#endif
                 }
                 break;

               case AllocTypeMatrixElement:
                 {
                   ElementPtr  pElement;
                   unsigned int i, tcount, flag;
                   tcount=ListPtr->AllocatedSize/sizeof(struct MatrixElement);
                   pElement=(ElementPtr)ListPtr->AllocatedPtr;
                   for (i=0; i<tcount; i++)
                   {
                     err|=MarkAlloc(Matrix,pElement->NextInRow,1);
                     err|=MarkAlloc(Matrix,pElement->NextInCol,1);
                     #if INITIALIZE
                       err|=MarkAlloc(Matrix,pElement->pInitInfo,1);
                     #endif
                     pElement++;
                   }
                 }
                 break;

               case AllocTypeFillinListNodeStruct:
               {
                 struct FillinListNodeStruct  *pElement;
                 ASSERT(ListPtr->AllocatedSize==sizeof(struct FillinListNodeStruct))
                 pElement=(struct FillinListNodeStruct  *)ListPtr->AllocatedPtr;
                 err|=MarkAlloc(Matrix,pElement->pFillinList,1);
                 err|=MarkAlloc(Matrix,pElement->Next,1);
               }
               break;

               case AllocTypeAllocationHashRecord:
               {
                 spAllocationListHashPtr  pElement;
                 unsigned int i,cnt;
//                 stat_AllocationRecord++;
                 pElement=(spAllocationListHashPtr)ListPtr->AllocatedPtr;
                 cnt=ListPtr->AllocatedSize/sizeof(struct spAllocationRecordHash);
                 for (i=0; i<cnt; i++)
                 {
                   unsigned int flag;
//                 err|=sysMarkAlloc(sys,pElement->AllocatedPtr,0);
                   err|=MarkAlloc(Matrix,pElement->NextRecord,1);
                   err|=MarkAlloc(Matrix,pElement->PrevRecord,1);
                   err|=MarkAlloc(Matrix,pElement->NextPeerRecord,1);
                   err|=MarkAlloc(Matrix,pElement->PrevPeerRecord,1);
                   err|=MarkAlloc(Matrix,pElement->AllocRecord,1);
                   pElement++;
                 }

            }
            break;

#if BUILDHASH

            case AllocTypespHtab:
            {
                struct spHtab *tab;
                unsigned int i, tcount;

                 tcount=ListPtr->AllocatedSize/sizeof(struct spHtab);
                 tab=(struct spHtab *)ListPtr->AllocatedPtr;

                 for (i=0; i<tcount; i++)
                 {
                   err|=MarkAlloc(Matrix,tab->entries,1);
                   tab++;
                 }
            }
            break;

            case AllocTypespHelt:
            {
                struct spHelt *elt;
                unsigned int i, tcount;

                 tcount=ListPtr->AllocatedSize/sizeof(struct spHelt);
                 elt=(struct spHelt *)ListPtr->AllocatedPtr;

                 for (i=0; i<tcount; i++)
                 {
                   err|=MarkAlloc(Matrix,elt->eptr,1);
                   err|=MarkAlloc(Matrix,elt->next,1);
                   elt++;
                 }
            }
            break;

            case AllocTypespHeblk:
            {
                struct heblk *eblk;
                unsigned int i, tcount;

                 tcount=ListPtr->AllocatedSize/sizeof(struct heblk);
                 eblk=(struct heblk *)ListPtr->AllocatedPtr;

                 for (i=0; i<tcount; i++)
                 {
                   unsigned int j;

                   for (j=0; j<HEBLKSZ; j++)
                   {
                      err|=MarkAlloc(Matrix,eblk->elts[j].eptr,1);
                      err|=MarkAlloc(Matrix,eblk->elts[j].next,1);

                    }
                   err|=MarkAlloc(Matrix,eblk->next,1);
                   eblk++;
                 }
            }
            break;

#endif // BUILDHASH

               case AllocTypeAllocationRecord:
               {
                 spAllocationListPtr  pElement;
                 unsigned int i, tcount;
                 tcount=ListPtr->AllocatedSize/sizeof(struct spAllocationRecord);
                 pElement=(spAllocationListPtr)ListPtr->AllocatedPtr;

                 for (i=0; i<tcount; i++)
                 {
                   err|=MarkAlloc(Matrix,pElement->AllocatedPtr,0);
                   err|=MarkAlloc(Matrix,pElement->NextRecord,1);
                   err|=MarkAlloc(Matrix,pElement->PrevRecord,1);
                   pElement++;
                 }
               }
               break;


            default:
              fprintf(stderr,"Unknown Allocation type %u\n",ListPtr->AllocatedType);
            }
            if  (err>0)
            {
              printf("Block allocated in %s: %u contains links outside of memory allocation table\n", ListPtr->AllocatedDebugFilename, ListPtr->AllocatedDebugLineNumber);
            }

        }
        ListPtr=ListPtr->NextRecord;
      }
    }
#if (1==1)
    /* Display orphant blocks and their size */
    allocs=0;
    allocsize=0;
    ListPtr=Matrix->TopOfAllocationList;
    while (ListPtr!=NULL)
    {
       if (ListPtr->Mark==0)
       {
        allocs++;
        allocsize+=ListPtr->AllocatedSize;
        printf("Orphant memory block allocated in %s:%u\n",ListPtr->AllocatedDebugFilename, ListPtr->AllocatedDebugLineNumber);
       }
       ListPtr=ListPtr->NextRecord;
    }
    if (allocs>0)
    {
      printf("Orphant memory blocks %u, Size: %lu kB\n",allocs, allocsize>>10);
      return 1;
    }
    else
    {
      printf ("No orphant found\n");
    }
    return 0;
#endif
}

#define spValidatePtr(Matrix,ptr) spValidatePtrFunc(Matrix,ptr,__FILE__,__LINE__)
void spValidatePtrFunc (MatrixPtr Matrix, void* ptr, const char* filename, unsigned int linenumber)
{
    spAllocationListPtr  ListPtr;
    if ((Matrix==NULL) || (ptr==NULL))
    {
       return;
    }
    ListPtr=Matrix->TopOfAllocationList;
    while (ListPtr!=NULL)
    {
        if (((unsigned long)ptr>=(unsigned long)ListPtr->AllocatedPtr) &&
            ((unsigned long)ptr<(((unsigned long)ListPtr->AllocatedPtr)+ListPtr->AllocatedSize)))
        {
            return;
        }
        ListPtr=ListPtr->NextRecord;
    }

    fprintf(stderr, "Pointer (0x%08lx) not found in allocation list %s:%u\n", (unsigned long)ptr, filename, linenumber);
    exit (-1);
    return;
}

