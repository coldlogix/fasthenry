/* This is the main part of the code */

#include "induct.h"
#include "sparse/spMatrix.h"
#include <string.h>
#include <ctype.h>
#include <sys/time.h>

#include "Prec_cost.h"
#include "BreakupSeg.h"
#include "SetupMulti.h"
#include "addgroundplane.h"
#include "barnoldi.h"
#include "calcp.h"
#include "cx_ludecomp.h"
#include "default_opts.h"
#include "deg_mutual.h"
#include "dist_betw_fils.h"
#include "findpaths.h"
#include "hole.h"
#include "joelself.h"
#include "mulSetup.h"
#include "fillM.h"
#include "gmres.h"
#include "SetupComputePsi.h"
#include "mutual.h"
#include "newPrecond.h"
#include "Precond.h"
#include "parse_command_line.h"
#include "readGeom.h"
#include "regurgitate.h"
#include "savemat_mod.h"
#include "writefastcap.h"
#include "memmgmt.h"
#include "gp.h"

#ifdef MPI
#include <mpi.h>
#endif

#define MAXPRINT 1
#define MAXCHARS 400
#define TIMESIZE 10

FILE *fp, *fp2, *fp3, *fb;
int num_exact_mutual;
int num_fourfil;
int num_mutualfil;
int num_found;
int num_perp;
int forced = 0;  /* for debugging inside exact_mutual() */

char outfname[200];
char outfname2[200];

/* savemat_mod machine type */
#ifdef DEC
int machine = 0000;
#else
int machine = 1000;
#endif

/* SRW */
charge *assignFil(SYS*, SEGMENT*, int*, charge*);
double **MatrixAlloc(int, int, int);
void fillA(SYS*);
void old_fillM(SYS*);
void fillZ(SYS*);
#if SUPERCON == ON
void fillZ_diag(SYS*, double);
void set_rvals(SYS*, double);
#endif
double resistance(FILAMENT*, double);
/* int matherr(struct exception*); */
int countlines(FILE*);
static int local_notblankline(char*);
void savemats(SYS*);
void savecmplx(FILE*, char*, CX**, int, int);
void savecmplx2(FILE*, char*, CX**, int, int);
void formMZMt(SYS*);
void oldformMZMt(SYS*);
void formMtrans(SYS*);
void compare_meshes(MELEMENT*, MELEMENT*);
void cx_dumpMat_totextfile(FILE*, CX**, int, int);
void dumpMat_totextfile(FILE*, double**, int, int);
void dumpVec_totextfile(FILE*, double*, int);
void fillMrow(MELEMENT**, int, double*);
void dump_to_Ycond(FILE*, int, SYS*);
void saveCarray(FILE*, char*, double**, int, int);
int nnz_inM(MELEMENT**, int);
void dump_M_to_text(FILE*, MELEMENT**, int, int);
void dump_M_to_matlab(FILE*, MELEMENT**, int, int, char*);
void pick_ground_nodes(SYS*);
int pick_subset(strlist*, SYS*);
void concat4(char*, char*, char*, char*);

/* wao: Refacturated functions from main */
void InitIndsys(SYS* indsys, int argc, char **argv);
void ReadGeometry(SYS* indsys,charge** chglist);
ssystem* MultipoleSetup(SYS* indsys, charge* chglist);
void Scanninggraph (SYS* indsys);
void OpenOutputfile (SYS* indsys, unsigned int MPIrank);
void AllocFunc1 (SYS* indsys);
void FillingM (SYS* indsys);
void CalcROM(SYS* indsys, ssystem* sys, charge* chglist);
void Precondition (int m,double freq, SYS* indsys,ssystem* sys);
void FormMtZM(int m, double freq, SYS* indsys);
void CalcConductor (int m, double freq, EXTERNAL* ext, int i, CX *b, CX *x0, CX* vect, SYS* indsys, ssystem* sys, charge* chglist, unsigned int MPIrank);


#define sysCommBcastMPI (0)
#define sysCommtoMPI (1)
#define sysCommfromMPI (2)
#define sysCommtoFile (3)
#define sysCommfromFile (4)

#define sysCopyMatrixoverMPI(ptr,rows,cols,size, grp) sysCopyMatrix(ptr,rows,cols,size,sysCommBcastMPI, (FILE*)&grp)
#define sysSendMatrixoverMPI(ptr,rows,cols,size,dest) sysCopyMatrix(ptr,rows,cols,size,sysCommtoMPI, (FILE*)&dest)
#define sysRecvMatrixoverMPI(ptr,rows,cols,size,src) sysCopyMatrix(ptr,rows,cols,size,sysCommfromMPI, (FILE*)&src)
#define sysMatrixSave(ptr,rows,cols,size, file) sysCopyMatrix(ptr,rows,cols,size,sysCommtoFile, file)
#define sysMatrixLoad(ptr,rows,cols,size, file) sysCopyMatrix(ptr,rows,cols,size,sysCommfromFile, file)
void sysCopyMatrix (double*** ptr, int rows, int cols, int size, int CommMode, FILE* fop);


#define sysCopyindsysoverMPI(ptr,grp) sysCopyindsys(ptr, sysCommBcastMPI, (FILE*)&grp)
#define sysSendindsysoverMPI(ptr,dest) sysCopyindsys(ptr, sysCommtoMPI, (FILE*)&dest)
#define sysRecvindsysoverMPI(ptr,src) sysCopyindsys(ptr, sysCommfromMPI, (FILE*)&src)
#define sysindsysSave(ptr,file) sysCopyindsys(ptr, sysCommtoFile, file)
#define sysindsysLoad(ptr,file) sysCopyindsys(ptr, sysCommfromFile, file)
void sysCopyindsys (SYS** psys, int CommMode, FILE* fop);

void sysDestroy( SYS** sys );
void sysDestroyMatrix (SYS* indsys, double*** ptr, int rows);
sysAllocationListHashPtr sysFindRecordAllocation (SYS* sys, char *Ptr );
void sysAnalyze2 (SYS* indsys);
unsigned int sysCalcHash(void* ptr);


typedef struct
{
   int CurrPreCond;
} MPIMachineState;
typedef MPIMachineState *MPIMachineStatePtr;

typedef struct _MPITask
{
  int m;
  int conductor;        /* -1: Calc preconditioner, 0...n Calculate conductor n) */
  double freq;
  int assigned;
  int PreCondSrc;      /* Source machine of preconditioner results, -1: It is already there */
  struct _MPITask *next;   /* Pointer to next structure */
  struct _MPITask **prev;  /* Pointer to next pointer in the previous structure */
} MPITask;
typedef MPITask *MPITaskPtr;

typedef struct {
  unsigned int m;
  unsigned int x;
  unsigned int y;
  CX data;
} MPIresult;


int main(int argc, char **argv)
{

  double width, height, length, freq, freqlast;
  int Linc, Winc, Hinc, filnum;
  double ratio;
  int i,j,k,m, last, err;
  char fname[80], tempstr[10];
  double r_height, r_width;
  CX dumb;
  CX *vect =0, **pvect =0;         /* space needed by gmres */
  double tol = 1e-8;
  double ftimes[TIMESIZE];
  CX *b =0, *x0 =0;
  SYS *indsys =0;           /* holds all the big variables (inductance system) */
  double totaltime;
  time_t t1,t2;
#ifdef MPI
  MPI_Comm MPI_WrkGrp;
#endif

  int MPIrank, MPIsize;
  int MPIdone;
  unsigned int MPImaxparallel;

  MPIMachineStatePtr MPIMachineStates;
  MPITaskPtr MPITasks;
  MPITaskPtr* MPICurrTask;
  unsigned int* MPICurrTransCnt;
  MPITaskPtr MPIsendbuf;
  MPITaskPtr MPIrecvbuf;
  CX*** MPIresults;

  const MPITask MPIdummyTask= {-1, -1, 0, 0, -1, NULL, NULL};

  t1=time(NULL);

  /* Initialize memory counters */
  memmgmt_init();

  /* Initialize MPI */
  MPIMachineStates=NULL;
  MPITasks=NULL;
  MPICurrTask=NULL;
  MPICurrTransCnt=NULL;

  MPIsize=1;
  MPIrank=0;
#ifdef MPI
  MPI_Init (&argc, &argv);	/* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);	/* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &MPIsize);	/* get number of processes */
#endif

  MPIdone=0;

  MPIsendbuf=NULL;
  if (MPIrank==0)
  {
    MPIsendbuf=(MPITaskPtr)malloc(sizeof(MPITask)*MPIsize);
  }
  MPIrecvbuf=(MPITaskPtr)malloc(sizeof(MPITask));


  for(i = 0; i < TIMESIZE; i++)
  {
    ftimes[i] = 0;
  }

  fb=NULL;

  if (MPIrank==0)
  {

  starttimer;
  indsys = (SYS *)calloc(1, sizeof(SYS));
  sysRecordAllocation(indsys,(char*)indsys, sizeof(SYS), sysAllocTypeindsys);

  sysALLOC(indsys->opts,1,ind_opts,ON,IND,indsys,sysAllocTypeindopts);

  InitIndsys(indsys, argc, argv);
  /* Register in Memory management */

  /* degugging counters for calls to functions used in mutual() */
  num_exact_mutual = 0;
  num_fourfil = 0;
  num_mutualfil = 0;
  num_found = 0;
  num_perp = 0;

  /** Read geometry **/
  ReadGeometry (indsys, &indsys->chglist);
  stoptimer;
  ftimes[0] = dtime;

   /** Multipole setup **/
  starttimer;
  MultipoleSetup(indsys, indsys->chglist);
  stoptimer;
  ftimes[5] = dtime;
  if (indsys->opts->debug == ON)
  {
    printf("Time for Multipole Setup: %lg\n",dtime);
  }

 /** Scanning graph **/
  starttimer;
  Scanninggraph (indsys);
  stoptimer;
  ftimes[6] = dtime;

  /** Open result file **/
  OpenOutputfile (indsys, MPIrank);

  /** Allocate memory **/
  AllocFunc1 (indsys);

  /** Form M and Z **/
  starttimer;
  FillingM (indsys);
  stoptimer;
  ftimes[1] = dtime;
  if (indsys->opts->debug == ON)
  {
    printf("Time to Form M and Z: %lg\n",dtime);
  }

  printf("Total Memory allocated: %ld kilobytes\n",memcount/1024);

  if (indsys->opts->debug == ON)
  {
    printf("Memory used and freed by lookup table: %d kilobytes\n",
	   get_table_mem());
  }

  /* free memory for lookup table */
  destroy_table();

  /*
     This function may alters
     indsys->precond_type
     indsys->precond_subtype
  */
  choose_and_setup_precond(indsys);

  if (indsys->opts->dumpMats)
  {
    printf("saving some files disk...\n");
    savemats(indsys);
  }

  if (indsys->logofstep == 0.0)
  {
    printf("no frequency range data read!\n");
    exit(1);
  }

  if (indsys->fmin == 0)
  {
      printf("***First frequency is zero. Only the DC case will be run.***\n");
      if (!indsys->dont_form_Z)
      {
        printf("Warning: First frequency is zero, but -sludecomp was not specified.\n\
      Use this setting to save time and memory.\n");
      }
  }

  if (indsys->opts->debug == ON)
  {
    /* open Ycond.mat */
    concat4(outfname,"Ycond",indsys->opts->suffix,".mat");
    /* SRW -- this is binary data */
    fp = fopen(outfname, "wb");
    if (fp == NULL)
    {
      printf("couldn't open file %s\n",outfname);
      exit(1);
    }
  }

  /* open b.mat */
  if (indsys->opts->dumpMats != OFF)
  {
    #ifdef MPI
      sprintf(outfname,"b_%u_%s.mat",MPIrank,indsys->opts->suffix);
    #else
      concat4(outfname,"b",indsys->opts->suffix,".mat");
    #endif
    /* SRW -- this is binary data */
    if ((fb = fopen(outfname,"wb")) == NULL)
    {
      fprintf(stderr, "No open fb\n");
    }
  }

   starttimer;
   CalcROM(indsys,indsys->sys,indsys->chglist);
   stoptimer;
   ftimes[8] += dtime;


  }

  /* Determine calculation tasks on root */
  if (MPIrank==0)
  {
    unsigned int i;
    EXTERNAL *ext;
    MPITaskPtr *ptptr;
    MPImaxparallel=0;

    MPIMachineStates=malloc(sizeof(MPIMachineState)*MPIsize);
    for (i=0; i<MPIsize;i++)
    {
      MPIMachineStates[i].CurrPreCond=-1;
    }

    ptptr=&MPITasks;
    for(freq = indsys->fmin, m = 0;
      (indsys->fmin != 0 && freq <= indsys->fmax*1.001) || (indsys->fmin == 0 && m == 0);
       m++, freq = (indsys->fmin != 0 ? pow(10.0,log10(indsys->fmin) + m*indsys->logofstep) : 0.0))
    {
      MPITaskPtr tptr;
      tptr=calloc(1,sizeof(MPITask));
      tptr->prev=ptptr;
      *ptptr=tptr;
      ptptr=&tptr->next;
      tptr->m=-1;
      tptr->freq=freq;
      tptr->conductor=m;
      tptr->assigned=-1;

      for(ext = get_next_ext(indsys->externals), i=0; ext != NULL;
                             ext = get_next_ext(ext->next),i++)
      {
        tptr=calloc(1,sizeof(MPITask));
        tptr->prev=ptptr;
        *ptptr=tptr;
        ptptr=&tptr->next;

        tptr->m=m;
        tptr->freq=freq;
        tptr->conductor=i;
        tptr->assigned=-1;
        MPImaxparallel++;

      }
    }
    MPICurrTask=calloc(MPIsize, sizeof(MPITaskPtr));
    MPICurrTransCnt=calloc(MPIsize, sizeof(unsigned int));
    MPIresults=calloc(m, sizeof(CX***));

    if (MPImaxparallel<MPIsize)
    {
      printf("WARNING: Utilize only %u of %u processors\n",MPImaxparallel,MPIsize);
    }
  }


#if 1==0
  sysAnalyze2 (indsys);
  printf("Rank %u check integrity indsys\n",MPIrank); fflush(stdout);
  sysAllocationAnalysis (indsys);
#endif


#ifdef MPI
  starttimer;
  MPI_Bcast(&MPImaxparallel, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* If there are more Nodes than possible parallel computations */
  /* Shrink communication group */
  {
    unsigned int WrkGrpMembers, *WrkGrpMembersList, i;
    MPI_Group GroupAll, MyGroup;

   WrkGrpMembers=MIN(MPImaxparallel,MPIsize);

   MPI_Comm_group(MPI_COMM_WORLD, &GroupAll);

   if (MPIrank>=WrkGrpMembers)
   {
     MPI_Group_incl(GroupAll, 1, &MPIrank, &MyGroup);
   }
   else
   {
     WrkGrpMembersList=malloc(WrkGrpMembers*sizeof(int));
     for (i=0; i<WrkGrpMembers; i++)
     {
       WrkGrpMembersList[i]=i;
     }

     MPI_Group_incl(GroupAll, WrkGrpMembers, WrkGrpMembersList, &MyGroup);
     free (WrkGrpMembersList);
     WrkGrpMembersList=NULL;
   }

   MPI_Comm_create(MPI_COMM_WORLD, MyGroup, &MPI_WrkGrp);

   MPI_Group_free(&MyGroup);

   /* Dismiss not used nodes */
   if (MPIrank>=WrkGrpMembers)
   {
     MPI_Finalize();
     exit (0);
   }
   MPIsize=WrkGrpMembers;

   /* Copy main data structure if we are more nodes than one */
   if (MPIsize>1)
   {
     printf("Rank %u Copy indsys over instances\n",MPIrank); fflush(stdout);
     sysCopyindsysoverMPI (&indsys, MPI_WrkGrp);
   }
   stoptimer;
   ftimes[7] = dtime;
  }
#endif

  {
    printf("Rank %u Allocate memory for solving\n",MPIrank); fflush(stdout);
    /* let's always dump the residual history */
    indsys->resids = sysMatrixAlloc(indsys, indsys->num_sub_extern,indsys->opts->maxiters, sizeof(double));
    indsys->resid_real = sysMatrixAlloc(indsys,indsys->num_sub_extern, indsys->opts->maxiters,
				     sizeof(double));
    indsys->resid_imag = sysMatrixAlloc(indsys,indsys->num_sub_extern, indsys->opts->maxiters,
				     sizeof(double));
    sysALLOC(indsys->niters,indsys->num_sub_extern,double,ON,IND,indsys,sysAllocTypeGeneric);
    indsys->FinalY = (CX **)sysMatrixAlloc(indsys,indsys->num_extern, indsys->num_sub_extern, sizeof(CX));

    sysALLOC(b    ,indsys->num_mesh,CX,ON,IND,indsys,sysAllocTypeGeneric);
    sysALLOC(x0,   indsys->num_mesh,CX,ON,IND,indsys,sysAllocTypeGeneric);
    sysALLOC(vect ,indsys->num_mesh,CX,ON,IND,indsys,sysAllocTypeGeneric);
    sysALLOC(pvect,indsys->num_mesh,CX*,ON,IND,indsys,sysAllocTypePtrArray);

    if (MPIrank>0)
    {
      /** Open result file **/
      OpenOutputfile (indsys, MPIrank);
    }
  }

  /* All machines enter main loop */
  while (MPIdone==0)
  {
    EXTERNAL *ext;
    int currm;

    /* plan the task assignment on root */
    if (MPIrank==0)
    {
      unsigned int i;

      /* Determine tasks for the machines */
      for (i=0; i<MPIsize; i++)
      {
        MPICurrTransCnt[i]=0;
        ASSERT(MPICurrTask[i]==NULL);
      }
      for (i=0; i<MPIsize; i++)
      {
        int found=0;
        MPITaskPtr ttask;

        /* Continue with existing frequency or do initial precond */
        ttask=MPITasks;
        found=0;
        while ((found==0) && (ttask!=NULL))
        {
          if ((ttask->m==MPIMachineStates[i].CurrPreCond) && (ttask->assigned==-1))
          {
            MPICurrTask[i]=ttask;
            ttask->assigned=i;
            ttask->PreCondSrc=-1;
            found=1;
          }
          ttask=ttask->next;
        }
      }
      for (i=0; i<MPIsize; i++)
      {
        int found=0;
        MPITaskPtr ttask;

        /* In case no task is assigned a) try to start with next frequency (cheap .. no additional data exchg required) */
        if ((MPICurrTask[i]==NULL) && (MPIMachineStates[i].CurrPreCond>-1))
        {
           // Reset preconditioner
           MPIMachineStates[i].CurrPreCond=-1;

          /* Continue with existing frequency or do initial precond */
          ttask=MPITasks;
          found=0;
          while ((found==0) && (ttask!=NULL))
          {
            if ((ttask->m==MPIMachineStates[i].CurrPreCond) && (ttask->assigned==-1))
            {
              MPICurrTask[i]=ttask;
              ttask->assigned=i;
              ttask->PreCondSrc=-1;
              found=1;
            }
            ttask=ttask->next;
          }
        }
      }
      for (i=0; i<MPIsize; i++)
      {
        int found=0;
        MPITaskPtr ttask;

        /* In case no task is assigned b) Copy PreConditioner results from other machine and take over conductor calculation */
        if (MPICurrTask[i]==NULL)
        {
          MPITaskPtr ttask_candidate;
          unsigned int ttask_transcnt;
          unsigned int ttask_candidate_src;

          ttask=MPITasks;
          ttask_candidate=NULL;
          while (ttask!=NULL)
          {
            if ((ttask->assigned==-1) && (ttask->m>-1))
            {
              /* This task is not assigned yet */
              unsigned int j;

              /* Has any machine the Preconditioner results ? */
              j=0;
              while (j<MPIsize)
              {
                if (MPIMachineStates[j].CurrPreCond==ttask->m)
                {
                  /* there is a candiate */
                  unsigned int takeit;
                  takeit=0;

                  /* Take first finding */
                  if (ttask_candidate==NULL)
                  {
                    takeit=1;
                  }
                  else
                  {
                    /* if following finding is better in terms of communication work load perfer this one over the other*/
                    if (MPICurrTransCnt[j]<ttask_transcnt)
                    {
                      takeit=1;
                    }
                  }

                  if (takeit>0)
                  {
                    ttask_candidate=ttask;
                    ttask_candidate_src=j;
                    ttask_transcnt=MPICurrTransCnt[ttask_candidate_src];
                  }
                }
                j++;
              }

            }
            ttask=ttask->next;
          }
          if (ttask_candidate!=NULL)
          {
            MPICurrTask[i]=ttask_candidate;
            ttask_candidate->PreCondSrc=ttask_candidate_src;
            ttask_candidate->assigned=i;
            MPICurrTransCnt[ttask_candidate_src]++;
          }

        }

      } /* for i=0...MPIsize */

      /* Are we done ? ==> No task assigned anymore */
      {
        unsigned int found;

        i=0;
        found=0;
        while ((i<MPIsize) && (found==0))
        {
          if (MPICurrTask[i]!=NULL)
          {
            found=1;
          }
          i++;
        }

        if (found==0)
        {
           /* All tasks shall be completed */
           ASSERT (MPITasks==NULL)
           MPIdone=1;
        }
      }
      /* Print game plan */
      {
        unsigned int i;
        printf("------");
        for (i=0; i<MPIsize; i++)
        {
          printf("------", MPIMachineStates[i].CurrPreCond);
        }
        printf("\n");
        printf("Cond: ");
        for (i=0; i<MPIsize; i++)
        {
          printf("%5i ", MPIMachineStates[i].CurrPreCond);
        }
        printf("\n");

        printf("Tran: ");
        for (i=0; i<MPIsize; i++)
        {
           if (MPICurrTask[i]!=NULL)
           {
              if (MPICurrTask[i]->PreCondSrc>-1)
              {
                printf("->%3i ", MPICurrTask[i]->PreCondSrc);
              }
              else
              {
                printf("      ");
              }
           }
           else
           {
             printf("      ");
           }
        }
        printf("\n");

        printf("Calc: ");
        for (i=0; i<MPIsize; i++)
        {
           if (MPICurrTask[i]!=NULL)
           {
            if (MPICurrTask[i]->m==-1)
            {
              printf("%2u:P  ", MPICurrTask[i]->conductor);
            }
            else
            {
              printf("%2u:%2u ", MPICurrTask[i]->m, MPICurrTask[i]->conductor);
            }
           }
           else
           {
             printf("      ");
           }
        }
        printf("\n");

      }
   }

   /* Bring everyone up to speed in terms of their tasks */
   starttimer;
   #ifdef MPI
   MPI_Bcast(&MPIdone, 1, MPI_INT, 0, MPI_WrkGrp);
   #endif

   if (MPIdone==0)
   {
     /* Assemble send buf */
     if (MPIrank==0)
     {
       unsigned int i;
       MPITaskPtr tptr;

       ASSERT (MPIsendbuf!=NULL)
       tptr=MPIsendbuf;
       for (i=0; i<MPIsize; i++)
       {
          if (MPICurrTask[i]!=NULL)
          {
            memcpy(tptr, MPICurrTask[i], sizeof(MPITask));
          }
          else
          {
            memcpy(tptr, &MPIdummyTask, sizeof(MPITask));
          }
          tptr++;
       }
     }
     ASSERT (MPIrecvbuf!=NULL)
     #ifdef MPI
     MPI_Scatter(MPIsendbuf,sizeof(MPITask)/sizeof(char),MPI_CHAR,
                 MPIrecvbuf,sizeof(MPITask)/sizeof(char),MPI_CHAR,
                 0, MPI_WrkGrp);
     #else
     memcpy (MPIrecvbuf, MPIsendbuf, sizeof(MPITask));
     #endif
   }

   #ifdef MPI
   /* Exchange preconditioner results as planed */

   if (MPIdone==0)
   {
     unsigned int Comm;
     unsigned int *Comms = NULL;
     unsigned int *CommList = NULL;
     unsigned int *CommsList = NULL;
     unsigned int Comms_all;
     unsigned GrpCommBufSize;
     unsigned int *GrpCommBuf = NULL;
     MPI_Group GroupAll, MyGroup;
     MPI_Comm MyComm;

     /* Sum outgoing transfers pers node */
     GrpCommBufSize=0;
     if (MPIrank==0)
     {
       unsigned int i;

       /* Sum outgoing communication per node */
       Comms=calloc(sizeof(int),MPIsize);
       Comms_all=0;
       for (i=0; i<MPIsize; i++)
       {
          if (MPICurrTask[i]!=NULL)
          {
              if (MPICurrTask[i]->PreCondSrc>-1)
              {
                 Comms[MPICurrTask[i]->PreCondSrc]++;
                 Comms_all++;
              }
          }
       }

       /* If there are 1:n communications MPI groups are used */
       /* Identify the cases */
       {
         for (i=0; i<MPIsize; i++)
         {
           if (Comms[i]>0)
           {
             GrpCommBufSize+=2+Comms[i];
           }
         }
       }
     }

     Comm=0;
     /* Distribute to nodes */
     MPI_Scatter(Comms, 1, MPI_INT,
                 &Comm, 1, MPI_INT,
                 0, MPI_WrkGrp);
     MPI_Bcast( &GrpCommBufSize, 1, MPI_INT, 0, MPI_WrkGrp);

     /* Assemble groups for 1:n communication */
     if (GrpCommBufSize>0)
     {
       GrpCommBuf=calloc(GrpCommBufSize,sizeof(int));
       if (MPIrank==0)
       {
         unsigned int *tbuf, tcnt ,i;
         tbuf=GrpCommBuf;
         tcnt=0;
         for (i=0; i<MPIsize; i++)
         {
           if (Comms[i]>0)
           {
             unsigned int j;
             ASSERT (tcnt<GrpCommBufSize)
             *(tbuf++)=1+Comms[i];
             *(tbuf++)=i;
             tcnt+=2;
             for (j=0; j<MPIsize; j++)
             {
               if (MPICurrTask[j]!=NULL)
               {
                 if (MPICurrTask[j]->PreCondSrc==i)
                 {
                   ASSERT (tcnt<GrpCommBufSize)
                   *(tbuf++)=j;
                   tcnt++;
                 }
               }
             }
           }
         }
       }
     }

     /* Broad cast groups and their members */
     MPI_Bcast( GrpCommBuf, GrpCommBufSize, MPI_INT, 0, MPI_WrkGrp);

     /* On root free allocated space for destination list */
     if (MPIrank==0)
     {
       free(Comms);
       Comms=NULL;
     }

     {
       unsigned int tcnt=0;
       unsigned int *tbuf = GrpCommBuf;
       unsigned int *myGroupBuf=NULL;
       while (tcnt<GrpCommBufSize)
       {
         unsigned int GrpMembers,j;
         unsigned int *tGroupBuf;
         tGroupBuf=tbuf;
         GrpMembers=*(tbuf++);
         ASSERT (GrpMembers>0)
         tcnt++;
         for (j=0; j<GrpMembers; j++)
         {
           if ( (*(tbuf++))==MPIrank )
           {
             /* This Rank is a member of this group */
             myGroupBuf=tGroupBuf;
           }
           tcnt++;
         }
       }
       ASSERT(tcnt==GrpCommBufSize)
       if (myGroupBuf!=NULL)
       {
/*
           unsigned int i;
           char txt[200],txt1[200];
           sprintf (txt,"Rank %u: Forming group with ", MPIrank);
           for (i=0; i<myGroupBuf[0]; i++)
           {
             sprintf(txt1," %u", myGroupBuf[1+i]);
             strcat(txt, txt1);
           }
*/
         MPI_Comm_group(MPI_WrkGrp, &GroupAll);
         MPI_Group_incl(GroupAll, myGroupBuf[0], &myGroupBuf[1], &MyGroup);
         MPI_Comm_create(MPI_WrkGrp, MyGroup, &MyComm);
 //        sprintf(txt1," Communicator %x\n",MyComm);strcat(txt, txt1);printf(txt); fflush (stdout);
       }

       if ((GrpCommBufSize>0) && (myGroupBuf==NULL))
       {
         /* All nodes of original group MPI_WrkGrp need to call MPI_Comm_create */
         /* Nodes that do not communicate create their very own ( and unused ) groups just containing themselfs */
         MPI_Comm_group(MPI_WrkGrp, &GroupAll);
         MPI_Group_incl(GroupAll, 1, &MPIrank, &MyGroup);
         MPI_Comm_create(MPI_WrkGrp, MyGroup, &MyComm);
 //        printf("Rank %u: Create group for myself Communicator %x\n",MPIrank,MyComm); fflush (stdout);
       }

       /* Process outgoing preconditioner result */
       if ((myGroupBuf!=NULL) && (Comm>0))
       {
          unsigned int MPIgrprank;
          /* This is sender, so Communication rank 0 is expected */
          MPI_Comm_rank (MyComm, &MPIgrprank);
          ASSERT (MPIgrprank==0)

          /* Broadcast preconditioner results */
//          printf("Rank %u: Broadcast Preconditioner results\n", MPIrank); fflush (stdout);
          sysCopyMatrixoverMPI((double***)&indsys->MtZM, indsys->num_mesh, indsys->num_mesh, sizeof(CX),MyComm);
 //         printf("Rank %u: Broadcast Preconditioner spMatrix results\n", MPIrank); fflush (stdout);
          spCopyMatrixoverMPI(&indsys->sparMatrix, MyComm);
       }
       else
       {
         ASSERT(Comm==0)
       }

       /* Check for incomming preconditioner result */
       if (MPIrecvbuf->m>-1)
       {
         if (MPIrecvbuf->PreCondSrc>-1)
         {
           if (indsys->MtZM!=NULL)
           {
             sysDestroyMatrix (indsys, (double***)&indsys->MtZM, indsys->num_mesh);
           }
           if (indsys->sparMatrix!=NULL)
           {
              spDestroy(indsys->sparMatrix);
              indsys->sparMatrix=NULL;
           }
           ASSERT (myGroupBuf!=NULL)
           {
             unsigned int MPIgrprank;

             /* This is receiver, so Group rank >0 is expected */
             /* Negative values indicate this node is not a member of the group - not expected*/
             MPI_Comm_rank (MyComm, &MPIgrprank);
             ASSERT (MPIgrprank>0)

//             printf("Rank %u: Receive Bcast Preconditioner results from %u\n", MPIrank, MPIrecvbuf->PreCondSrc); fflush(stdout);
             sysCopyMatrixoverMPI((double***)&indsys->MtZM, indsys->num_mesh, indsys->num_mesh, sizeof(CX),MyComm);
//             printf("Rank %u: Receive Bcast Preconditioner spMatrix results from %u\n", MPIrank, MPIrecvbuf->PreCondSrc); fflush(stdout);
             spCopyMatrixoverMPI(&indsys->sparMatrix, MyComm);

           }
         }
       }

       if (GrpCommBufSize>0)
       {
//         printf("Rank %u: Release Comm and Group\n", MPIrank);
         MPI_Comm_free(&MyComm);
         MPI_Group_free(&MyGroup);
       }
     }

     if (GrpCommBuf!=NULL)
     {
       free(GrpCommBuf);
       GrpCommBuf=NULL;
     }


   }

   if (MPIrank==0)
   {
     unsigned int i;
     for (i=0; i<MPIsize; i++)
     {
       if (MPICurrTask[i]!=NULL)
       {
         if (MPICurrTask[i]->PreCondSrc>-1)
         {
           MPIMachineStates[i].CurrPreCond=MPIMachineStates[MPICurrTask[i]->PreCondSrc].CurrPreCond;
         }
       }
     }
   }
   /* Exchange preconditioner results as planed */
   #endif
   stoptimer;
   ftimes[7] += dtime;

   /* Execute assigned task */
   ext=NULL;
   currm=-1;
   if (MPIdone==0)
   {
     if (MPIrecvbuf->m==-1)
     {
        if (MPIrecvbuf->conductor>-1)
        {
          m=MPIrecvbuf->m;
          freq=MPIrecvbuf->freq;
          printf("Frequency = %lg\n",freq); fflush(stdout);

          starttimer;
          FormMtZM(m, freq, indsys);
          stoptimer;
          ftimes[2] += dtime;

          starttimer;
          Precondition (m, freq, indsys, indsys->sys);
          stoptimer;
          ftimes[3] += dtime;

          if (indsys->opts->debug == ON)
          {
            printf("Rank %u: Time spent on forming Precond: %lg\n",MPIrank, dtime);
          }

          //spAnalyzeAllocatedMemory(indsys->sparMatrix);

       }
     }
     else
     {

       unsigned int i;

       currm=m=MPIrecvbuf->m;
       freq=MPIrecvbuf->freq;
       for(ext = get_next_ext(indsys->externals), i=0;
            (ext != NULL) && (i!=MPIrecvbuf->conductor);
            ext = get_next_ext(ext->next),i++);

       starttimer;
       CalcConductor(m, freq,
                     ext, i,
                     b, x0, vect,
                     indsys,indsys->sys, indsys->chglist, MPIrank );
       stoptimer;
       ftimes[4] += dtime;

     }

     if (MPIrank==0)
     {
        unsigned int i;
        for (i=0; i<MPIsize; i++)
        {
           if (MPICurrTask[i]!=NULL)
           {
             if (MPICurrTask[i]->m==-1)
             {
               MPIMachineStates[i].CurrPreCond=MPICurrTask[i]->conductor;
             }
           }
        }
     }
   }
   /* End of execute assigned task */

   /* Send computation result back to root */
   starttimer;
   {
     int elems, elems_all;
     int* MPIrecvelems = NULL;
     int* tlens = NULL;
     int* toffs = NULL;
     MPIresult* tresult = NULL;
     MPIresult* tresults= NULL;

     /* Count of results on each node */
     elems=0;
     if (ext!=NULL)
     {
       EXTERNAL *text;

       for(text = indsys->externals; text != NULL; text = text->next)
       {
         elems++;
       }
     }

     /* Allocate memory on root node for storage of count of results of each node */
     if (MPIrank==0)
     {
       MPIrecvelems=malloc(sizeof(int)*MPIsize);
     }
     /* Gather count of results of each node on root node */
     #ifdef MPI
     MPI_Gather(&elems, 1, MPI_INT,
                MPIrecvelems, 1, MPI_INT,
                0, MPI_WrkGrp);
     #else
      MPIrecvelems[0]=elems;
     #endif
     /* Create of result buffer of this node to be sent to root node */
     if (elems>0)
     {
       EXTERNAL *text;
       MPIresult* tptr;
       ASSERT(currm>-1);
       tresult=malloc(sizeof(MPIresult)*elems);
       tptr=tresult;
       for(text = indsys->externals; text != NULL; text = text->next)
       {
         tptr->m=currm;
         tptr->x=text->Yindex;
         tptr->y=ext->col_Yindex;
         memcpy(&tptr->data,&indsys->FinalY[tptr->x][tptr->y], sizeof(CX));
         tptr++;
       }
     }

     /* Allocate memory on root node for all results */
     if (MPIrank==0)
     {
        int tcnt;

        /* Calculate length and offset of response of all nodes in result buffer */
        tlens = malloc(sizeof(int)*MPIsize);
        toffs = malloc(sizeof(int)*MPIsize);
        tcnt=0;
        elems_all=0;
        for (i=0; i<MPIsize; i++)
        {

          toffs[i]=tcnt;
          tlens[i]=MPIrecvelems[i]*sizeof(MPIresult)/sizeof(char);
          tcnt+=tlens[i];

          elems_all+=MPIrecvelems[i];

        }
        /* Allocate memory for results of all nodes */
        tresults=calloc(elems_all,sizeof(MPIresult));

     }

     /* Gather results from all nodes */
     #ifdef MPI
     MPI_Gatherv(tresult,  (sizeof(MPIresult)*elems)/sizeof(char), MPI_CHAR,
                 tresults, tlens, toffs, MPI_CHAR,
                 0, MPI_WrkGrp);
     #else
     memcpy(tresults, tresult, elems*sizeof(MPIresult));
     #endif

     /* Process result on root node */
     if (MPIrank==0)
     {
       MPIresult* tptr;

       /* Free length and offsets in receive buffer */
       free (tlens); tlens=NULL;
       free (toffs); toffs= NULL;

       /* Process results from all nodes */
       tptr=tresults;
       for (i=0;i<elems_all; i++)
       {
   //      printf("Result: %u %u %u\n",tptr->m, tptr->x,tptr->y); fflush (stdout);

         /* Create result matrix if not there yet */
         if (MPIresults[tptr->m]==0)
         {
           MPIresults[tptr->m]=(CX **)sysMatrixAlloc(indsys,indsys->num_extern, indsys->num_sub_extern, sizeof(CX));
         }
         /* store data in final data matrix */
         memcpy(&MPIresults[tptr->m][tptr->x][tptr->y],&tptr->data, sizeof(CX));

         /* Process next element */
         tptr++;
       }

       /* Free transmission buffer */
       free (tresults);
       tresults=NULL;
     }
     stoptimer;
     ftimes[7] += dtime;


     /* Release memory for results of each node */
     if (elems>0)
     {
       free (tresult);
       tresult=NULL;
     }
   }
   /* End of send computation result back to root */


   /* Remove completed tasks from chain */
   if (MPIrank==0)
   {
        unsigned int i;
        for (i=0; i<MPIsize; i++)
        {

           if (MPICurrTask[i]!=NULL)
           {
             /* Remove from chain */
             if (MPICurrTask[i]->next!=NULL)
             {
               MPICurrTask[i]->next->prev=MPICurrTask[i]->prev;
             }
             *(MPICurrTask[i]->prev)=MPICurrTask[i]->next;

             /* Free memory */
             free(MPICurrTask[i]);
             MPICurrTask[i]=NULL;
           }

        }
//        printf("\n");
    }

  } /* while MPIdone */

  if (MPIrank>0)
  {
    free (MPIsendbuf);
    MPIsendbuf=NULL;
    free(MPICurrTask);
    MPICurrTask=NULL;
    free(MPICurrTransCnt);
    MPICurrTransCnt=NULL;
  }
  free (MPIrecvbuf);
  MPIrecvbuf=NULL;

  if (MPIrank==0)
  {
      for(freq = indsys->fmin, m = 0;
          (indsys->fmin != 0 && freq <= indsys->fmax*1.001) || (indsys->fmin == 0 && m == 0);
           m++, freq = (indsys->fmin != 0 ? pow(10.0,log10(indsys->fmin) + m*indsys->logofstep) : 0.0))
      {
           if (indsys->opts->debug == ON)
           {
              /* dump matrix to matlab file */
              static char fname[20], tempstr[5];

              sprintf(tempstr, "%d", m);

              strcpy(fname,"Ycond");
              strcat(fname,tempstr);
              savecmplx(fp, fname, MPIresults[m], indsys->num_extern,indsys->num_sub_extern);
           }

           /* dump matrix to text file */
           if (indsys->num_extern == indsys->num_sub_extern)
           {
                fprintf(fp3, "Impedance matrix for frequency = %lg %d x %d\n ", freq,
	            indsys->num_extern, indsys->num_extern);
                cx_invert(MPIresults[m], indsys->num_extern);
           }
           else
           {
                fprintf(fp3, "ADMITTANCE matrix for frequency = %lg %d x %d\n ", freq,
	            indsys->num_extern, indsys->num_sub_extern);
           }

           cx_dumpMat_totextfile(fp3, MPIresults[m], indsys->num_extern, indsys->num_sub_extern);
           fflush(fp3);
      }
   }

  if (fp!=NULL)
  {
    fclose(fp);
  }

  if (fp3!=NULL)
  {
    fclose(fp3);
  }

  if (fb!=NULL)
  {
    printf("Rank %u Close DumpMats\n",MPIrank); fflush(stdout);
    fclose(fb);
  }



  if (MPIrank==0)
  {
    printf("\nAll impedance matrices dumped to file Zc%s.mat\n\n",indsys->opts->suffix);

    if (indsys->opts->debug == ON)
    {
      printf("Calls to exact_mutual: %15d\n",num_exact_mutual);
      printf("         fourfils:     %15d\n",num_fourfil);
      printf("         mutualfil:    %15d\n",num_mutualfil);
      printf("Number found in table: %15d\n",num_found);
      printf("Number perpendicular:  %15d\n",num_perp);
      printf("\n");
    }
  }
  {
      double gtimes[TIMESIZE];
      unsigned int i;

      #ifdef MPI
      MPI_Reduce(ftimes, gtimes, 10, MPI_DOUBLE, MPI_SUM, 0,MPI_WrkGrp );
      #else
      for (i=0; i<TIMESIZE; i++)
      {
        gtimes[i]=ftimes[i];
      }
      #endif
      totaltime = 0;
      for(i = 0; i < TIMESIZE; i++)
      {
        totaltime += gtimes[i];
      }

      if (MPIrank==0)
      {
         double elapsedTime;

         printf("Times:  Read geometry   %lg\n",gtimes[0]);
         printf("        Multipole setup %lg\n",gtimes[5]);
         printf("        Scanning graph  %lg\n",gtimes[6]);
         printf("        Form A M and Z  %lg\n",gtimes[1]);
         printf("        form M'ZM       %lg\n",gtimes[2]);
         printf("        Form precond    %lg\n",gtimes[3]);
         printf("        GMRES time      %lg\n",gtimes[4]);
         printf("        MPI             %lg\n",gtimes[7]);
         printf("        ROM             %lg\n",gtimes[8]);
         printf("   Total:               %lg\n",totaltime);

         t2=time(NULL);
         elapsedTime = t2-t1;
         printf("   Total real time      %lg\n",elapsedTime);

      }

      #ifdef MATTDEBUG
        /* print memory bins */
        for(i = 0; i < 1001; i++)
        {
          fprintf(stderr, "%d\n", membins[i]);
        }
      #endif
  }

  sysDestroy(&indsys);


#ifdef MPI
   /* All ranks */
   MPI_Finalize();

#endif

  return (0);
}

void InitIndsys(SYS* indsys, int argc, char **argv)
{
  /* Init global structure */

  indsys->externals = NULL;
  indsys->num_extern = 0;
  indsys->segment = NULL;
  indsys->nodes = NULL;
  indsys->planes = NULL;                               /* CMS 7/2/92 */
  indsys->logofstep = 0;
  Parse_Command_Line(indsys, indsys->opts, argc, argv);

  /*indsys->r_height = indsys->r_width = indsys->opts->filratio;*/

  /* set the type equal to that specified on command line for now */
  indsys->precond_type = indsys->opts->precond;
}

void ReadGeometry (SYS* indsys,charge** chglist)
{
  FILE* fp;
  int err;
  char outfname[200];
  char outfname2[200];

  if (indsys->opts->fname == NULL)
  {
    printf("No input file given.  Reading from stdin...\n");
    fp = stdin;
  }
  else
  {
    if (strcmp(indsys->opts->fname,"-") == 0)
    {
      printf("Reading from stdin...\n");
      fp = stdin;
    }
    else
    {
      /* SRW -- read ascii file */
      fp = fopen(indsys->opts->fname, "r");
      if (fp == NULL)
      {
        printf("Couldn't open %s\n", indsys->opts->fname);
        exit(1);
      }
      printf("Reading from file: %s\n",indsys->opts->fname);
    }
  }
  /* read in geometry */
  err = readGeom(fp, indsys);
  if (err != 0)
  {
    exit(err);
  }
  fclose(fp);

  /*  fprintf(stdout,"Number of nodes before breaking up: %d\n",
      indsys->num_nodes); */

  /* regurgitate input file to stdout */
  if (indsys->opts->regurgitate == TRUE)
  {
    regurgitate(indsys);
  }

  if ((indsys->opts->makeFastCapFile & HIERARCHY)
      && !(indsys->opts->makeFastCapFile & (SIMPLE | REFINED)))
  {
    /* the hierarchy has been dumped already, and nothing more to do. so quit*/
    exit(0);
  }

  /* make a file suitable for keith's postscript maker */
  if (indsys->opts->makeFastCapFile & SIMPLE)
  {
    fprintf(stdout, "Making simple zbuf file...");
    fflush(stdout);
    concat4(outfname,"zbuffile",indsys->opts->suffix,"");
    concat4(outfname2,outfname,"_","shadings");
    writefastcap(outfname, outfname2, indsys);
    fprintf(stdout, "Done\n");
    if (! (indsys->opts->makeFastCapFile & REFINED) )
    {
      /* we are done, after making fastcap files for visualization */
      exit(0);
    }
  }

  /* initialize Gaussian quadrature arrays for mutual(). */
  /* See dist_betw_fils.c */
  fill_Gquad();

  /* initialize Freeable allocation stuff for mutual terms lookup table */
  init_table();

  /* break each segment into filaments */
  {
    charge *chgend, chgdummy;
    SEGMENT *seg;

    indsys->num_fils = 0;
    chgend = &chgdummy;
    for(seg = indsys->segment; seg != NULL; seg = seg->next)
    {
      chgend = assignFil(indsys, seg, &indsys->num_fils, chgend);
    }
    chgend->next = NULL;
    *chglist = chgdummy.next;
  }

}

ssystem* MultipoleSetup(SYS* indsys, charge* chglist)
{
  ssystem* sys;

#if SUPERCON == ON
  set_rvals(indsys, 2*PI*indsys->fmin);
#endif

  /* set up multipole stuff */
  sys = SetupMulti(chglist, indsys);

#if 1 == 0
  if (indsys->opts->debug == ON && 1 == 0)  /* for debugging eval pass */
    dump_evalcnts(sys);
#endif
 return sys;
}

void Scanninggraph (SYS* indsys)
{
  printf("Scanning graph to find fundamental circuits...\n");
  /* find all the necessary meshes */
  make_trees(indsys);

  /* clear node->examined and seg->is_deleted */
  clear_marks(indsys);

  /* find meshes in groundplane due to holes and put them at the front
     of the list of trees so that they are created by fillM() first so
     that they are marked before regular ground plane meshes  */
  find_hole_meshes(indsys);

  /* determine if a subset of columns is to be computed and assign
     col_Yindex.  This will happen if -x option is used. */
  indsys->num_sub_extern = pick_subset(indsys->opts->portlist, indsys);
}

void OpenOutputfile (SYS* indsys, unsigned int MPIrank)
{
   char outfname[200];
  EXTERNAL *ext;

  /* Write to Zc.mat only if we aren't only running for visualization */
  if ( !(indsys->opts->makeFastCapFile & (SIMPLE | REFINED))
       && !(indsys->opts->orderROM > 0 && indsys->opts->onlyROM)   )
  {
    if (MPIrank>0)
      sprintf(outfname,"Zc%s_%u.mat",indsys->opts->suffix, MPIrank);
    else
      concat4(outfname,"Zc",indsys->opts->suffix,".mat");   /* put filnames together */
    /* SRW -- this is ascii data */
    fp3 = fopen(outfname, "w");
    if (fp3 == NULL)
    {
      printf("couldn't open file %s\n",outfname);
      exit(1);
    }

    for(ext = indsys->externals; ext != NULL; ext=ext->next)
    {
      /* printf("Row %d :  %s  to  %s\n",ext->Yindex,ext->source->node[0]->name,
         ext->source->node[1]->name); */

      if (ext->portname[0] == '\0')
      {
        fprintf(fp3,"Row %d:  %s  to  %s\n",ext->Yindex+1,ext->name1,ext->name2);
      }
      else
      {
        fprintf(fp3,"Row %d:  %s  to  %s, port name: %s\n",
                ext->Yindex+1,ext->name1,ext->name2,ext->portname);
      }
    }

    if (indsys->num_sub_extern != indsys->num_extern)
    {
      for(ext = indsys->externals; ext != NULL; ext=ext->next)
      {
        if (ext->col_Yindex != -1)
        {
          fprintf(fp3,"Col %d: port name: %s\n",
                  ext->col_Yindex+1,ext->portname);
        }
      }
    }
  }

}

void AllocFunc1 (SYS* indsys)
{

/*  FILE *fp3; */

  int last,i;
  int nonp, planemeshes;           /* CMS 6/7/92 */
  GROUNDPLANE *plane;                   /* CMS 7/7/92 */
  NODES *node;
  EXTERNAL *ext;


  /* count nodes */
  indsys->num_real_nodes = 0;
  indsys->num_nodes = 0;
  last = -1;
  for(node = indsys->nodes, i = 0; node != NULL; node = node->next, i++)
  {
    indsys->num_nodes++;
    if (getrealnode(node) == node)
    {
      indsys->num_real_nodes++;
      if (last + 1 != i)
      {
      /* printf("Non equivalent nodes must be listed first??,
	 no I take that back.\n"); */
      /* exit(1); */
      }
      last = i;
    }
  }

  /* CMS 7/3/92 ------------------------------------------------------------*/
  planemeshes = 0;
  nonp = 0;

  for(plane = indsys->planes; plane != NULL; plane = plane->next)
  {
    planemeshes = planemeshes + plane->numesh;
    if(plane->external == 0)
    {
      nonp++;
    }
  }
  /*------------------------------------------------------------------------*/

/* moved lower for shading
  if (indsys->opts->makeFastCapFile & REFINED) {
    fprintf(stdout, "Making refined zbuf file...");
    fflush(stdout);

    concat4(outfname,"zbuffile2",indsys->opts->suffix,"");


    concat4(outfname2,outfname,"_","shadings");
    writefastcap(outfname, outfname2, indsys);
    fprintf(stdout, "Done\n");
  }
*/

  /* figure out number of meshes */
  indsys->tree_meshes = count_tree_meshes(indsys->trees);

  indsys->extra_meshes = indsys->tree_meshes;

#if 1==0
  unimplemented junk
  indsys->extra_meshes = estimate_extra_meshes(indsys->trees, FILS_PER_MESH);
#endif
                                                          /* CMS 7/2/92 */
  indsys->num_mesh = indsys->num_fils - indsys->num_segs + indsys->extra_meshes + planemeshes;;
  if (count_externals(indsys->externals) != indsys->num_extern)
  {
    fprintf(stderr, "main:  discrepancy in num_extern and actual number\n");
    exit(1);
  }

  /*  indsys->dont_form_Z = indsys->fmin == 0 && indsys->opts->mat_vect_prod == DIRECT
                        && indsys->opts->soln_technique == ITERATIVE; */

  indsys->dont_form_Z = indsys->fmin == 0 && indsys->opts->soln_technique == LUDECOMP;


  printf("Number of Groundplanes : %d \n",indsys->num_planes);        /* CMS 7/2/92 */
  printf("Number of filaments: %10d\n", indsys->num_fils);
  printf("Number of segments:  %10d\n", indsys->num_segs);
  printf("Number of nodes:     %10d\n", indsys->num_nodes);
  printf("Number of meshes:    %10d\n", indsys->num_mesh);
  printf("          ----from tree:                   %d \n",indsys->tree_meshes);
  printf("          ----from planes: (before holes)  %d \n",planemeshes);
                                                            /* CMS 7/6/92 */
  printf("Number of conductors:%10d   (rows of matrix in Zc.mat) \n",
	 indsys->num_extern);
  printf("Number of columns:   %10d   (columns of matrix in Zc.mat) \n",
	 indsys->num_sub_extern);

  printf("Number of real nodes:%10d\n", indsys->num_real_nodes);  /* non equivalent */

/*  A = MatrixAlloc(indsys->num_real_nodes, indsys->num_fils, sizeof(double)); */

  sysALLOC(indsys->meshsect,indsys->num_extern+nonp+1,int,ON,IND,indsys,sysAllocTypeGeneric);
  sysALLOC(indsys->m_info,indsys->num_mesh,Minfo,ON,IND,indsys,sysAllocTypeGeneric);

  sysALLOC(indsys->Precond,indsys->num_mesh,PRE_ELEMENT*,ON,IND,indsys,sysAllocTypePtrArray);

  sysALLOC(indsys->Mlist,indsys->num_mesh,MELEMENT*,ON,IND,indsys,sysAllocTypePtrArray);
  sysALLOC(indsys->Mtrans,indsys->num_fils,MELEMENT*,ON,IND,indsys,sysAllocTypePtrArray);

  sysALLOC(indsys->R,indsys->num_fils,double,ON,IND,indsys,sysAllocTypeGeneric);

  if (indsys->opts->mat_vect_prod == DIRECT || indsys->opts->soln_technique == LUDECOMP)
  {
    if (!indsys->dont_form_Z)
    {
      indsys->Z = sysMatrixAlloc(indsys, indsys->num_fils, indsys->num_fils, sizeof(double));
    }
  }

/*
  if (indsys->opts->mat_vect_prod == DIRECT || indsys->opts->soln_technique == LUDECOMP)
  {
    for(i = 0; i < indsys->num_fils; i++)
    {
      for(j = 0; j < indsys->num_mesh; j++)
      {
	    indsys->M[i][j] = 0;
      }
    }
  }
*/
}

void FillingM (SYS* indsys)
{
  char outfname[200];
  char outfname2[200];
  int num_mesh;
  int MPIrank;

/*
  printf("filling A...\n");
  fillA(segment, A, indsys->num_segs);
*/

  printf("filling M...\n");
  num_mesh=indsys->num_mesh;
  fillM(indsys);  /* expects hole meshes to be marked */
                /* This function updates indsys->num_mesh */

  if (indsys->opts->makeFastCapFile & REFINED)
  {
    fprintf(stdout, "Making refined zbuf file...");
    fflush(stdout);
    concat4(outfname,"zbuffile2",indsys->opts->suffix,"");

    /* shadings file */
    concat4(outfname2,outfname,"_","shadings");

    /* output file to later turn into .ps. Note: this affects seg->is_deleted*/
    writefastcap(outfname, outfname2, indsys);
    fprintf(stdout, "Done\n");
    exit(0);  /* we are done after making fastcap files for visualization */;
  }


  if (num_mesh != indsys->num_mesh)
  {
    printf("Exact number of meshes after ground plane holes removed: %d\n",
	   indsys->num_mesh);
  }
  /* let's save a little space and allocate MZMt now. */
  if (!indsys->dont_form_Z
      && (indsys->opts->mat_vect_prod == DIRECT || indsys->opts->soln_technique == LUDECOMP))
  {
    indsys->MtZM = (CX **)sysMatrixAlloc(indsys, indsys->num_mesh, indsys->num_mesh, sizeof(CX));
    printf("Rank %u Alloc MtZM %08lx %u %u\n", MPIrank,(unsigned long)indsys->MtZM ,indsys->num_mesh, sizeof(CX)); fflush (stdout);
  }

  if (indsys->opts->dumpMats & MESHES)
  {
    if (indsys->opts->kind & MATLAB)
    {
      dump_mesh_coords(indsys);
    }
    else
    {
      dump_ascii_mesh_coords(indsys);
    }
  }

  formMtrans(indsys); /* Form M transpose by row*/

  printf("filling R and L...\n");
  fillZ(indsys);

}

void CalcROM(SYS* indsys, ssystem* sys, charge* chglist)
{
  /* Model Order Reduction: create, compute and print the model if requested */
  if (indsys->opts->orderROM > 0)
  {
    double *dtemp;

    FILE *fROM;
    double **B, **C;
    char *MRMt;            /* may be needed if ROM is requested */
    int actual_order;
    EXTERNAL *ext;
    int j;

    actual_order=0;

    /* create the matrix whose inverse we will need */
    createMRMt(&MRMt, indsys);
    /* now create the input and output matrices for the state */
    B = sysMatrixAlloc(indsys,indsys->num_mesh, indsys->num_extern, sizeof(double));
    C = sysMatrixAlloc(indsys,indsys->num_mesh, indsys->num_extern, sizeof(double));
    for(ext = indsys->externals, j = 0;
        ext != NULL; ext = ext->next, j++)
    {
      int_list *elem;
      /* create C first; C = N, where Vs = N Vterminal */
      for(elem = ext->indices; elem != NULL; elem = elem->next)
      {
        if (elem->index > indsys->num_mesh || ext->Yindex > indsys->num_extern)
        {
	       fprintf(stderr, "Indexing into matrix C out of bounds\n");
	       exit(1);
	    }
        C[elem->index][ext->Yindex] = elem->sign;
        B[elem->index][ext->Yindex] = elem->sign;
      }
    }

    if (indsys->opts->mat_vect_prod != MULTIPOLE && !indsys->dont_form_Z)
    {
      printf("multiplying M*(L)*transpose(M)for model order reduction\n");
      /* this form of storage isn't the best */
      formMLMt(indsys);      /*form (M^t)*(L)*M and store in indys->MtZM*/


      actual_order = ArnoldiROM(B, C, (double **)NULL, (char**)MRMt, indsys->num_mesh,
                                indsys->num_extern, indsys->num_extern, indsys->opts->orderROM,
                                realMatVect, indsys, sys, chglist);
    }
    else
    {
      if (indsys->opts->mat_vect_prod == MULTIPOLE)
      {
        actual_order = ArnoldiROM(B, C, (double **)NULL, (char**)MRMt, indsys->num_mesh,
                                indsys->num_extern, indsys->num_extern, indsys->opts->orderROM,
                                realComputePsi, indsys, sys, chglist);
      }
    }

    if (indsys->opts->debug == ON)
    {
#if 1==0
      /* open orig.mat */
      concat4(outfname,"orig",indsys->opts->suffix,".mat");
      /* SRW -- this is binary data */
      if ((fROM = fopen(outfname,"wb")) == NULL)
        fprintf(stderr, "No open fROM\n");
      /* dump what we have of the original system */
      dumpROMbin(fROM, NULL, B, C, NULL,
                 indsys->num_mesh, indsys->num_extern, indsys->num_extern);
      fclose(fROM);
#endif

    /* open rom.m */
      concat4(outfname,"rom",indsys->opts->suffix,".m");
      /* SRW -- this is ascii data */
      if ((fROM = fopen(outfname,"w")) == NULL)
      {
        fprintf(stderr, "No open fROM\n");
      }

    /* now dump the reduced order model */
      dumpROM(fROM, indsys->Ar, indsys->Br, indsys->Cr, indsys->Dr,
              actual_order * indsys->num_extern, indsys->num_extern, indsys->num_extern);
      /* close and save away the file */
      fclose(fROM);
    }

    /* generate equivalent circuit */
    concat4(outfname,"equiv_circuitROM",indsys->opts->suffix,".spice");
    /* SRW -- this is ascii data */
    if ((fROM = fopen(outfname,"w")) == NULL)
    {
      fprintf(stderr, "No open fROM\n");
    }
    /* now dump the reduced order model */
    dumpROMequiv_circuit(fROM, indsys->Ar, indsys->Br, indsys->Cr, indsys->Dr,
            actual_order * indsys->num_extern, indsys->num_extern, indsys->num_extern,
            indsys->title, indsys->opts->suffix, indsys);
    /* close and save away the file */
    fclose(fROM);
    /* end of Model Order Reduction generation */
  }

  /* NOTE: exit here if all we wnat is the reduced order model (ROM) */
  if (indsys->opts->orderROM > 0 && indsys->opts->onlyROM)
  {
    exit(0);
  }
}

/* uses sys->directlist,sys->length, sys->minx, sys->miny, sys->minz */
/*       indsys->Mtrans, indsys->Mlist, indsys->R, indsys->Z  */
/* Generates indsys->Precond */
/* Generates indsys->sparMatrix */
void Precondition (int m, double freq,  SYS* indsys,ssystem* sys)
{
  int err;
  char outfname[200];

  switch (indsys->opts->soln_technique)
  {
    case ITERATIVE:
      switch (indsys->precond_type)
      {
        case LOC:  /* ITERATIVE, LOC */
        {
	      printf("Forming local inversion preconditioner\n");

          switch (indsys->opts->mat_vect_prod)
          {
            case DIRECT: /* ITERATIVE, LOC, DIRECT */
              /* uses sys->directlist */
              /*       indsys->Mtrans, indsys->Mlist */
              /* Generates indsys->Precond */
              indPrecond_direct(sys, indsys, 2*PI*freq);
              break;

            case MULTIPOLE: /* ITERATIVE, LOC, MULTIPOLE */
              /* uses sys->directlist,sys->length, sys->minx, sys->miny, sys->minz */
              /*       indsys->Mtrans, indsys->Mlist, indsys->R, indsys->Z  */
              /* Generates indsys->Precond */

	          indPrecond(sys, indsys, 2*PI*freq);
	          break;
	        default:
	          fprintf(stderr, "Internal error: mat_vect_prod == %d\n",
		             indsys->opts->mat_vect_prod);
	          exit(1);
          }
        }
        break;

        case SPARSE: /* SPARSE */
        {
	      printf("Forming sparse matrix preconditioner..\n");
          /* uses sys->directlist,sys->length, sys->minx, sys->miny, sys->minz, sys->cubes */
          /*       indsys->Mtrans, indsys->Mlist, indsys->R, indsys->Z, indsys->segment  */
          /* Generates indsys->Precond */
          /* Generates indsys->sparMatrix */
          create_sparMatrix(indsys);
   	      fill_spPre(sys, indsys, 2*PI*freq);

          if (m == 0)
          {
	        if (indsys->opts->debug == ON)
	        {
	           printf("Number of nonzeros before factoring: %d\n",
		            spElementCount(indsys->sparMatrix));
            }

            /* Reorder and Factor the matrix */
	        err = spOrderAndFactor(indsys->sparMatrix, NULL, 1e-3, 0.0, 1);

	        if (indsys->opts->debug == ON)
	        {
              printf("Number of fill-ins after factoring: %d\n",
		          spFillinCount(indsys->sparMatrix));
            }
	      }
	      else
	      {
	        err = spFactor(indsys->sparMatrix);
	      }


	      if (err != 0)
	      {
	        fprintf(stderr,"Error on factor: %d\n",err);
	        exit(1);
	      }
        }
        break;

        default:
	        fprintf(stderr,"Unknown precond type: %d\n",indsys->precond_type);
	        exit(1);

      }
      break;

    case LUDECOMP: /* LUDECOMP */

      if (!indsys->dont_form_Z)
      {
        printf("Performing LU Decomposition...\n");
        cx_ludecomp(indsys->MtZM, indsys->num_mesh, FALSE);
        printf("Done.\n");
      }
      else
      {

        /* let's form a sparse version since fmin=0 and the L matrix = 0 */
        /* Generates indsys->sparMatrix */
        create_sparMatrix(indsys);
        printf("Filling sparse version M R Mt...");

        /* Alters indsys->sparMatrix */
        /* uses  indsys->Mtrans,indsys->R, indsys->segment  */
        fill_diagR(indsys);

        /* dump matrix to disk if requested */
        if (indsys->opts->dumpMats & MZMt)
        {
          printf("saving sparse MZMt to MZMt.mat...\n");
          concat4(outfname,"MZMt",indsys->opts->suffix,".mat");
          if (spFileMatrix(indsys->sparMatrix, outfname, "MZMt", 0, 1, 1) == 0)
          {
            fprintf(stderr,"saving sparse matrix failed\n");
          }
        }

        if (indsys->opts->debug == ON)
        {
          printf("Number of nonzeros before factoring: %d\n",
                 spElementCount(indsys->sparMatrix));
        }
        /* Reorder and Factor the matrix */
        /* since w = 0, this matrix is real but all complex computations
           are done since that is how the sparse package was compiled.
           But changing it doesn't help that much.  Must be dominated by
           reordering and such.*/
        printf("Factoring the sparse matrix...\n");
        err = spOrderAndFactor(indsys->sparMatrix, NULL, 1e-3, 0.0, 1);

        printf("Done factoring.\n");
        if (indsys->opts->debug == ON)
        {
          printf("Number of fill-ins after factoring: %d\n",
                 spFillinCount(indsys->sparMatrix));
        }
      }

      break;

    default:
        fprintf(stderr, "Internal error: Unknown type of soln_technique: %d\n",
	        indsys->opts->soln_technique);
        exit(1);
    }


}

void FormMtZM(int m, double freq, SYS* indsys)
{
  FILE *fp2;
  char outfname[200];

#if SUPERCON == ON
    if (freq != indsys->fmin)
      fillZ_diag(indsys, 2*PI*freq);
#endif

    if (!indsys->dont_form_Z
        && (indsys->opts->mat_vect_prod == DIRECT || indsys->opts->soln_technique==LUDECOMP))
    {
      printf("multiplying M*(R + jL)*transpose(M)\n");
      formMZMt(indsys);      /*form transpose(M)*(R+jL)*M  (no w) */

      if (indsys->opts->dumpMats & MZMt)
      {
	    if (m == 0)
	    {
	      if (indsys->opts->kind & MATLAB)
	      {
            concat4(outfname,"MZMt",indsys->opts->suffix,".mat");
            /* SRW -- this is binary data */
	        if ( (fp2 = fopen(outfname,"wb")) == NULL)
	        {
	          printf("Couldn't open file\n");
	          exit(1);
	        }
	        printf("Saving MZMt...\n");
	        savecmplx2(fp2,"MZMt",indsys->MtZM, indsys->num_mesh,indsys->num_mesh);
	        fclose(fp2);
	       }
	       if (indsys->opts->kind & TEXT)
	       {
             concat4(outfname,"MZMt",indsys->opts->suffix,".dat");
             /* SRW -- this is ascii data */
             fp2 = fopen(outfname,"w");
	         if ( fp2 == NULL)
	         {
	           printf("Couldn't open file\n");
	           exit(1);
	         }
	         cx_dumpMat_totextfile(fp2, indsys->MtZM,
				  indsys->num_mesh,indsys->num_mesh );
	         fclose(fp2);
	      }
	    }
      }

      printf("putting in frequency \n");

      /* put in frequency */
      {
        CX **MtZM;
        int i,j;

        MtZM = indsys->MtZM;
        for(i = 0; i < indsys->num_mesh; i++)
        {
	       for(j = 0; j < indsys->num_mesh; j++)
	       {
	          MtZM[i][j].imag *= 2*PI*freq;
           }
        }
      }
    }
}

void CalcConductor (int m, double freq,
                    EXTERNAL* ext, int i,
                    CX *b, CX *x0, CX* vect,
                    SYS* indsys, ssystem* sys, charge* chglist, unsigned int MPIrank)
{
    int j;
    int maxiters;
    char fname[80], tempstr[10];
     printf("Rank %u: conductor %d from node %s\n",MPIrank, i, get_a_name(ext->source)); fflush(stdout);

      /* initialize b */
      for(j = 0; j < indsys->num_mesh; j++)
      {
	    b[j] = x0[j] = CXZERO;
      }

      fill_b(ext, b);

      sprintf(fname, "b%d_%d",m,i);
      if (indsys->opts->dumpMats != OFF)
      {
	    savecmplx(fb, fname, &b, 1, indsys->num_mesh);
      }

      maxiters = MIN(indsys->opts->maxiters, indsys->num_mesh+2);

#if OPCNT == ON
      maxiters = 2;
#endif

      if (indsys->opts->soln_technique == ITERATIVE)
      {
	     printf("Calling gmres...\n");
	     if (indsys->opts->mat_vect_prod == MULTIPOLE)
	     {
//	        printf("%u SetupComputePsi \n",MPIrank);
	        gmres(indsys->MtZM, b, x0, inner, SetupComputePsi, indsys->num_mesh, maxiters, indsys->opts->tol,
		          sys, chglist, 2*PI*freq, indsys->R, indsys, i);
	     }
	     else
	     {
//            printf("%u DirectMatVec \n",MPIrank);
	        gmres(indsys->MtZM, b, x0, inner, directmatvec,    indsys->num_mesh, maxiters, indsys->opts->tol,
		          sys, chglist, 2*PI*freq, indsys->R, indsys, i);
	     }

	     if (indsys->precond_type == LOC)
         {
            int k;
	        multPrecond(indsys->Precond, x0, vect, indsys->num_mesh);
	        for (k=0; k < indsys->num_mesh; k++)
            {
	          x0[k] = vect[k];
            }
	     }
	     else
	     {
	       if (indsys->precond_type == SPARSE)
	       {
	         spSolve(indsys->sparMatrix, (spREAL*)x0, (spREAL*)x0);
	       }
	     }
      }
      else
      {
        if (!indsys->dont_form_Z)
        {
          cx_lu_solve(indsys->MtZM, x0, b, indsys->num_mesh);
        }
        else
        {
          spSolve(indsys->sparMatrix, (spREAL*)b, (spREAL*)x0);
        }
      }

      if (indsys->opts->dumpMats & GRIDS)
      {
	    /* Do stuff to look at groundplane current distribution */
	    makegrids(indsys, x0, ext->Yindex, m);
      }

      /* Extract appropriate elements of x0 that correspond to final Yc */
      extractYcol(indsys->FinalY, x0, ext, indsys->externals);

      if (indsys->opts->debug == ON)
      {
        FILE* fptemp;
        /* SRW -- this is binary data */
        #ifdef MPI
	      sprintf(tempstr, "Ytemp_%d", MPIrank);
	      fptemp = fopen(tempstr, "wb");
	    #else
	      fptemp = fopen("Ytemp.mat", "wb");
	    #endif
	    if (fptemp == NULL)
        {
	      printf("couldn't open file %s\n","Ytemp.mat");
	       exit(1);
	    }
	    strcpy(fname,"Ycond");
	    sprintf(tempstr, "%d", m);
        strcat(fname,tempstr);
	    savecmplx(fptemp, fname, indsys->FinalY, indsys->num_extern, indsys->num_sub_extern);
	    fclose(fptemp);
      }

}


/* this will divide a rectangular segment into many filaments */
charge *assignFil(SYS* indsys, SEGMENT *seg, int *num_fils, charge *chgptr)
{

  int i,j,k;
  double x,y,z, delw, delh;
  double hx, hy, hz;
  int temp, counter;
  int Hinc, Winc;
  double Hdiv, Wdiv;
  double height, width;
  double h_from_edge, w_from_edge, min_height, min_width;
  FILAMENT *tempfil, *filptr;
  int countfils;
  double wx, wy, wz, mag;  /* direction of width */
  NODES *node0, *node1;
  charge *clast, *chg;
  charge *ctemp;  /* temporary place so can allocate all the charges at once*/
  surface *dummysurf;
  int indices[4], row, col;
  double r_width, r_height;
            /* ratio of element sizes (for geometrically increasing size)*/

  Hinc = seg->hinc;
  Winc = seg->winc;
  r_height = seg->r_height;
  r_width = seg->r_width;

  clast = chgptr;

  seg->num_fils = Winc*Hinc;
  sysALLOC(seg->filaments,Hinc*Winc,FILAMENT,ON,IND,indsys,sysAllocTypeFilament);

  ctemp=NULL;
  sysALLOC(ctemp,Hinc*Winc,charge,ON,IND,indsys,sysAllocTypeCharge);

  dummysurf=NULL;
  sysALLOC(dummysurf,1,surface,ON,IND,indsys,sysAllocTypeSurface);
  dummysurf->type = CONDTR;

  /*  To make the filaments have geometrically decreasing areas */
  /* Hdiv and Wdiv are 1/(smallest Width) */

  if (fabs(1.0 - r_height) < EPS)
  {
    Hdiv = Hinc;
  }
  else
  {
    temp = Hinc/2;
    Hdiv = 2*(1.0 - pow(r_height, (double)temp))/(1.0 - r_height);
    if (Hinc%2 == 1) Hdiv += pow(r_height, (double)(temp));
  }

  if (fabs(1.0 - r_width) < EPS)
  {
    Wdiv = Winc;
  }
  else
  {
    temp = Winc/2;
    Wdiv = 2*(1.0 - pow(r_width, (double)temp))/(1.0 - r_width);
    if (Winc%2 == 1) Wdiv += pow(r_width, (double)(temp));
  }

  node0 = seg->node[0];
  node1 = seg->node[1];
  /* determine direction of width */
  if (seg->widthdir != NULL)
  {
    wx = seg->widthdir[XX];
    wy = seg->widthdir[YY];
    wz = seg->widthdir[ZZ];
  }
  else
  {
    /* default for width direction is in x-y plane perpendic to length*/
    /* so do cross product with unit z*/
    wx = -(node1->y - node0->y)*1.0;
    wy = (node1->x - node0->x)*1.0;
    wz = 0;
    if ( fabs(wx/seg->length) < EPS && fabs(wy/seg->length) < EPS)
    {
      /* if all of x-y is perpendic to length, then choose x direction */
      wx = 1.0;
      wy = 0;
    }
    mag = sqrt(wx*wx + wy*wy + wz*wz);
    wx = wx/mag;
    wy = wy/mag;
    wz = wz/mag;
  }
  /* height direction perpendicular to both */
  hx = -wy*(node1->z - node0->z) + (node1->y - node0->y)*wz;
  hy = -wz*(node1->x - node0->x) + (node1->z - node0->z)*wx;
  hz = -wx*(node1->y - node0->y) + (node1->x - node0->x)*wy;
  mag = sqrt(hx*hx + hy*hy + hz*hz);
  hx = hx/mag;
  hy = hy/mag;
  hz = hz/mag;

  filptr = seg->filaments;
  counter = 0;

  /* this will fill the 'filament' array. It uses symmetry wrt z and y */
  /* it generates the four corner fils first, then the next one in...  */
  /* 6/25/93 - added stuff to place fils in filament array so that adjacent */
  /*           fils are adjacent in the array.  This will make the meshes   */
  /*           in M consist of fils that are near each other.               */
  h_from_edge = 0.0;
  min_height = 1.0/Hdiv;  /* fil of smallest height */
  for(i = 0; i < Hinc/2 || (Hinc%2 == 1 && i == Hinc/2); i++)
  {

    /* height of the filaments for this row */
    height = (seg->height/Hdiv)*pow(r_height, (double)i);
    /* delh = (seg->height)*(0.5 - 1.0/Hdiv*
			  ( (1 - pow(r_height, (double)(i+1)))/(1-r_height)
			   - 0.5*pow(r_height, (double) i) )
			 );
			 */
    if (i == 0)
    {
      h_from_edge += min_height/2;
    }
    else
    {
      h_from_edge += min_height/2*pow(r_height,(double)(i-1))*(1+r_height);
    }
    delh = (seg->height)*(0.5 - h_from_edge);

    if (delh < 0.0 && fabs(delh/seg->height) > EPS)
    {
      printf("uh oh, delh < 0. delh/height = %lg\n", delh/seg->height);
    }

    w_from_edge = 0;
    min_width = 1.0/Wdiv;
    for(j = 0; j < Winc/2 || (Winc%2 == 1 && j == Winc/2); j++)
    {
      width = (seg->width/Wdiv)*pow(r_width, (double)j );
      /*delw = (seg->width)*
	(0.5 - 1.0/Wdiv*
	 ( (1 - pow(r_width, (double)(j+1)))/(1 - r_width)
	  - 0.5*pow(r_width, (double) j) )
	 );
	 */

      if (j == 0)
      {
	    w_from_edge += min_width/2;
	  }
      else
      {
	    w_from_edge += min_width/2*pow(r_width,(double)(j-1))*(1+r_width);
	  }

      delw = (seg->width)*(0.5 - w_from_edge);

      if (delw < 0.0 && fabs(delw/seg->width) > EPS)
      {
	    printf("uh oh, delw < 0. delw/width = %lg\n", delw/seg->width);
	  }
/*    tempfil = filptr; */
      countfils = 0;

      row = i;
      col = j;
      if (row%2 == 1)
      {
	    col = (Winc - 1) - col;
	  }
      indices[countfils] = col + Winc*row;
      filptr = &(seg->filaments[indices[countfils]]);

      filptr->x[0] = node0->x + hx*delh + wx*delw;
      filptr->x[1] = node1->x + hx*delh + wx*delw;
      filptr->y[0] = node0->y + hy*delh + wy*delw;
      filptr->y[1] = node1->y + hy*delh + wy*delw;
      filptr->z[0] = node0->z + hz*delh + wz*delw;
      filptr->z[1] = node1->z + hz*delh + wz*delw;
      filptr->pchg = &ctemp[counter++];
/*    filptr++; */
      countfils++;

      /* do symmetric element wrt y */
      if(j != Winc/2)
      {
	    row = i;
	    col = (Winc - 1) - j;
	    if (row%2 == 1)
	    {
	      col = (Winc - 1) - col;
	    }
	    indices[countfils] = col + Winc*row;
	    filptr = &(seg->filaments[indices[countfils]]);

	    filptr->x[0] = node0->x + hx*delh - wx*delw;
	    filptr->x[1] = node1->x + hx*delh - wx*delw;
	    filptr->y[0] = node0->y + hy*delh - wy*delw;
	    filptr->y[1] = node1->y + hy*delh - wy*delw;
	    filptr->z[0] = node0->z + hz*delh - wz*delw;
	    filptr->z[1] = node1->z + hz*delh - wz*delw;
	    filptr->pchg = &ctemp[counter++];
/*	filptr++; */
	    countfils++;
      }

      /* wrt z */
      if(i != Hinc/2)
      {
        row = (Hinc - 1) - i;
	    col = j;
	    if (row%2 == 1)
	    {
	      col = (Winc - 1) - col;
	    }
	    indices[countfils] = col + Winc*row;
	    filptr = &(seg->filaments[indices[countfils]]);

	    filptr->x[0] = node0->x - hx*delh + wx*delw;
	    filptr->x[1] = node1->x - hx*delh + wx*delw;
	    filptr->y[0] = node0->y - hy*delh + wy*delw;
	    filptr->y[1] = node1->y - hy*delh + wy*delw;
	    filptr->z[0] = node0->z - hz*delh + wz*delw;
	    filptr->z[1] = node1->z - hz*delh + wz*delw;
	    filptr->pchg = &ctemp[counter++];
	    filptr++;
	    countfils++;
      }

      /* wrt z and y */
      if( i != Hinc/2 && j != Winc/2)
      {
	    row = (Hinc - 1) - i;
	    col = (Winc - 1) - j;
	    if (row%2 == 1)
	      col = (Winc - 1) - col;
	    indices[countfils] = col + Winc*row;
	    filptr = &(seg->filaments[indices[countfils]]);

	    filptr->x[0] = node0->x - hx*delh - wx*delw;
	    filptr->x[1] = node1->x - hx*delh - wx*delw;
	    filptr->y[0] = node0->y - hy*delh - wy*delw;
	    filptr->y[1] = node1->y - hy*delh - wy*delw;
	    filptr->z[0] = node0->z - hz*delh - wz*delw;
	    filptr->z[1] = node1->z - hz*delh - wz*delw;
	    filptr->pchg = &ctemp[counter++];
/*	filptr++; */
	    countfils++;
      }

      for(k = 0; k < countfils; k++)
      {
	    tempfil = &(seg->filaments[indices[k]]);
	    tempfil->length = seg->length;
	    tempfil->area = width*height;
	    tempfil->width = width;
	    tempfil->height = height;
	    tempfil->filnumber = (*num_fils)++;
	    tempfil->segm = seg;

	    tempfil->lenvect[XX] = tempfil->x[1] - tempfil->x[0];
	    tempfil->lenvect[YY] = tempfil->y[1] - tempfil->y[0];
	    tempfil->lenvect[ZZ] = tempfil->z[1] - tempfil->z[0];
	    /* do stuff for multipole */
	    /*   make linked list entry */
	    chg = tempfil->pchg;
	    clast->next = chg;
	    clast = chg;
	    /* fill charge structure */
	    chg->max_diag = chg->min_diag = tempfil->length;
	    chg->x = (tempfil->x[0] + tempfil->x[1])/2.0;
	    chg->y = (tempfil->y[0] + tempfil->y[1])/2.0;
	    chg->z = (tempfil->z[0] + tempfil->z[1])/2.0;
	    chg->surf = dummysurf;
	    chg->dummy = FALSE;
	    chg->fil = tempfil;
/*	tempfil++; */
      }
    }
  }

  i = 0;
  while(i < Hinc*Winc && seg->filaments[i].pchg != NULL)
  {
    i++;
  }

  if (i != Hinc*Winc)
  {
    fprintf(stderr, "Hey, not all filaments created in assignfil()! \n");
    exit(1);
  }

  return clast;
}

double **MatrixAlloc(int rows, int cols, int size)
{

  double **temp;
  int i;

  temp = (double **)MattAlloc(rows,sizeof(double *));
  if (temp == NULL)
  {
    printf("not enough space for matrix allocation\n");
    exit(1);
  }

  for(i = 0; i < rows; i++)
  {
    temp[i] = (double *)MattAlloc(cols,size);
  }

  if (temp[rows - 1] == NULL)
  {
    printf("not enough space for matrix allocation\n");
    exit(1);
  }
  return(temp);
}

void fillA(SYS *indsys)
{
  SEGMENT *seg;
  NODES *node1, *node2, *node;
  MELEMENT **Alist;
  int i, counter;
  FILAMENT *fil;

  sysALLOC(indsys->Alist,indsys->num_real_nodes,MELEMENT*,ON,IND,indsys,sysAllocTypePtrArray);

  pick_ground_nodes(indsys);

  Alist = indsys->Alist;
  counter = 1;  /* ground is chosen already */

  for(seg = indsys->segment; seg != NULL; seg = seg->next)
  {
    node1 = getrealnode(seg->node[0]);
    node2 = getrealnode(seg->node[1]);
    if (node1->index == -1)
    {
      node1->index = counter++;
      Alist[node1->index] = NULL;
    }
    if (node2->index == -1)
    {
      node2->index = counter++;
      Alist[node2->index] = NULL;
    }
    for(i = 0; i < seg->num_fils; i++)
    {
      fil = &seg->filaments[i];
      Alist[node1->index] = insert_in_list(indsys,make_melement(indsys,fil->filnumber,
							 fil, 1),
					   Alist[node1->index]);
      Alist[node2->index] = insert_in_list(indsys,make_melement(indsys, fil->filnumber,
							 fil, -1),
					   Alist[node2->index]);
    }
  }

  if (counter != indsys->num_real_nodes - indsys->num_trees + 1)
  {
    fprintf(stderr,"Internal error when forming A: counter %d != num_real_nodes %d\n",
	    counter, indsys->num_real_nodes);
  }

}

/* this fills the kircoff's voltage law matrix (Mesh matrix) */
/* it maps a matrix of mesh currents to branch currents */
/* it might actually be what some think of as the transpose of M */
/* Here, M*Im = Ib  where Im are the mesh currents, and Ib the branch */
/* 6/92 I added Mlist which is a vector of linked lists to meshes.
   This replaces M.  But I keep M around for checking things in matlab. */
void old_fillM(SYS *indsys)
{
}

void fillZ(SYS *indsys)
{
  int i, j, k, m;
  FILAMENT *fil_j, *fil_m;
  int filnum_j, filnum_m;
  double w;
  SEGMENT *seg1, *seg2;
  double **Z, *R, freq;
  int num_segs;

  Z = indsys->Z;
  R = indsys->R;

  for(seg1 = indsys->segment; seg1 != NULL; seg1 = seg1->next) {
    for(j = 0; j < seg1->num_fils; j++) {
      fil_j = &(seg1->filaments[j]);
      filnum_j = fil_j->filnumber;
#if SUPERCON == ON
      R[filnum_j] = fil_j->length*seg1->r1/fil_j->area;
#else
      R[filnum_j] = resistance(fil_j, seg1->sigma);
#endif
      if (indsys->opts->mat_vect_prod != MULTIPOLE
          && !indsys->dont_form_Z) {
	for(seg2 = indsys->segment; seg2 != NULL; seg2 = seg2->next) {
	  for(m = 0; m < seg2->num_fils; m++) {
	    fil_m = &(seg2->filaments[m]);
	    filnum_m = fil_m->filnumber;
	    if (filnum_m == filnum_j) {
	      Z[filnum_m][filnum_m] = selfterm(fil_m); /* do self-inductance */
#if SUPERCON == ON
	      if (seg1->lambda != 0.0)
	        Z[filnum_m][filnum_m] += seg1->r2*fil_j->length/fil_j->area;
#endif
	    }
	    else
	      if (filnum_m > filnum_j) /*we haven't done it yet */

		Z[filnum_m][filnum_j]
		  = Z[filnum_j][filnum_m] = mutual(fil_j, fil_m);
	  }
	}
      }  /* end if (multipole or dont_form_Z) */
    }
  }
}

#if SUPERCON == ON
/* Have to reset the diag elements for each omega, if superconductor
 * (lambda != 0) and sigma != 0.
 */
void fillZ_diag(SYS *indsys, double omega)
{
  int i, j;
  FILAMENT *fil_j;
  int filnum_j;
  SEGMENT *seg1;
  double **Z, *R;
  double tmp, tmp1, dnom, r1, r2;

  Z = indsys->Z;
  R = indsys->R;

  for(seg1 = indsys->segment; seg1 != NULL; seg1 = seg1->next) {
    if (seg1->lambda != 0.0 && seg1->sigma != 0.0) {
      /* segment is a superconductor with frequency dependent terms */
      tmp = MU0*seg1->lambda*seg1->lambda;
      tmp1 = tmp*omega;
      dnom = tmp1*seg1->sigma;
      dnom = dnom*dnom + 1.0;
      seg1->r1 = r1 = seg1->sigma*tmp1*tmp1/dnom;
      seg1->r2 = r2 = tmp/dnom;
      for(j = 0; j < seg1->num_fils; j++) {
        fil_j = &(seg1->filaments[j]);
        filnum_j = fil_j->filnumber;
        R[filnum_j] = fil_j->length*r1/fil_j->area;
        Z[filnum_j][filnum_j] = selfterm(fil_j) +
          r2*fil_j->length/fil_j->area;
      }
    }
  }
}

/*
 * Theory:
 *    In normal metal:     (1)  del X H = -i*omega*mu*sigma * H
 *    In superconductor:   (2)  del X H = (1/lambda)^2 * H
 *
 *    In fasthenry, (1) is solved, so the game is to replace sigma in (1)
 *    with a complex variable that includes and reduces to (2).  We choose
 *
 *    sigma_prime = sigma + i/(omega*mu*lambda^2)
 *
 *    Then, using sigma_prime in (1) rather than sigma, one obtains an
 *    expression that reduces to (2) at omega = 0,  yet retains properties
 *    of (1).  This is the two-fluid model, where the sigma in sigma_prime
 *    represents the conductivity due to unpaired electrons.
 *
 *    Since sigma_prime blows up at omega = 0, we work with the impedance,
 *    which we take as z = r1 + i*omega*r2 = i/sigma_prime.  The r1 and
 *    r2 variables are thus
 *
 *           (3) r1 =    sigma*(omega*mu*lambda^2)^2
 *                    --------------------------------
 *                    (sigma*omega*mu*lambda^2)^2 + 1
 *
 *           (4) r2 =           mu*lambda^2
 *                    --------------------------------
 *                    (sigma*omega*mu*lambda^2)^2 + 1
 */

/* Set the r1, r2 entries to the appropriate values.  The r1 entry
 * comes from the real value of sigma (1/sigma for normal metals).
 * The r2 entry comes from the imaginary part of sigma, which reduces
 * to MU0*lambda^2 when the real part of sigma (default 0 for
 * superconductors) is negligible, and is 0 for normal metals.
 * The r2 term contributes to the inductance matrix diagonal.
 */
void set_rvals(SYS *indsys, double omega)
{
  SEGMENT *seg1;
  double tmp, tmp1, dnom, r1, r2;

  for(seg1 = indsys->segment; seg1 != NULL; seg1 = seg1->next) {
    if (seg1->lambda != 0.0) {
      /* segment is a superconductor */
      tmp = MU0*seg1->lambda*seg1->lambda;
      if (seg1->sigma != 0.0) {
        /* the terms become frequency dependent */
        tmp1 = tmp*omega;
        dnom = tmp1*seg1->sigma;
        dnom = dnom*dnom + 1.0;
        r1 = seg1->sigma*tmp1*tmp1/dnom;
        r2 = tmp/dnom;
      }
      else {
        r1 = 0.0;
        r2 = tmp;
      }
    }
    else {
      r1 = 1/seg1->sigma;
      r2 = 0.0;
    }
    seg1->r1 = r1;
    seg1->r2 = r2;
  }
}
#endif

/* calculates resistance of filament */
double resistance(FILAMENT *fil, double sigma)
/* double sigma;  conductivitiy */
{
  return  fil->length/(sigma * fil->area);
}

/* mutual inductance functions moved to mutual.c */

/*
int matherr(struct exception *exc)
{
  printf("Err in math\n");
  return(0);
}
*/

/* This counts the nonblank lines of the file  fp (unused) */
int countlines(FILE *fp)
{
  int count;
  char temp[MAXCHARS], *returnstring;

  count = 0;
  while( fgets(temp,MAXCHARS, fp) != NULL)
    if ( local_notblankline(temp) ) count++;

  return count;
}

/* returns 1 if string contains a nonspace character */
static int local_notblankline(char *string)
{
   while( *string!='\0' && isspace(*string))
     string++;

   if (*string == '\0') return 0;
     else return 1;
}

/* This saves various matrices to files and optionally calls fillA() if
   the incidence matrix, A, is requested */
void savemats(SYS *indsys)
{
  int i, j;
  FILE *fp, *fp2;
  FILAMENT *tmpf;
  SEGMENT *tmps;
  double *dumb, *dumbx, *dumby, *dumbz;
  int num_mesh, num_fils, num_real_nodes;
  double *Mrow;   /* one row of M */
  double **Z, *R;
  int machine, counter, nnz, nnz0;
  ind_opts *opts;
  MELEMENT *mesh;
  MELEMENT **Mlist = indsys->Mlist;
  MELEMENT **Alist;

  num_mesh = indsys->num_mesh;
  num_fils = indsys->num_fils;
  num_real_nodes = indsys->num_real_nodes;

  Mrow = (double *)MattAlloc(num_fils,sizeof(double));
  R = indsys->R;
  Z = indsys->Z;
  opts = indsys->opts;

  if (opts->dumpMats & DUMP_M) {
    /* count non-zero entries */
    nnz = nnz_inM(indsys->Mlist, num_mesh);

    if (opts->kind & TEXT) {
      concat4(outfname,"M",opts->suffix,".dat");
      /* SRW -- this is ascii data */
      if ( (fp2 = fopen(outfname,"w")) == NULL) {
	printf("Couldn't open file\n");
	exit(1);
      }
      dump_M_to_text(fp2, Mlist, num_mesh, nnz);
      fclose(fp2);
    }

    if (opts->kind & MATLAB) {
      concat4(outfname,"M",opts->suffix,".mat");
      /* SRW -- this is binary data */
      if ( (fp = fopen(outfname,"wb")) == NULL) {
	printf("Couldn't open file\n");
	exit(1);
      }
      dump_M_to_matlab(fp, Mlist, num_mesh, nnz, "M");
      fclose(fp);
    }
  }

  if (opts->dumpMats & DUMP_A) {

    /* fill the A matrix */
    fillA(indsys);
    Alist = indsys->Alist;

    /* count non-zero entries */
    nnz = nnz_inM(indsys->Alist, num_real_nodes);

    /* how many non-zeros in ground node 0? */
    nnz0 = nnz_inM(indsys->Alist, 1);

    /* now dump it to a file without ground node */

    if (opts->kind & TEXT) {
      concat4(outfname,"A",opts->suffix,".dat");
      /* SRW -- this is ascii data */
      if ( (fp2 = fopen(outfname,"w")) == NULL) {
	printf("Couldn't open file\n");
	exit(1);
      }
      dump_M_to_text(fp2, &Alist[1], num_real_nodes - 1, nnz - nnz0);
      fclose(fp2);
    }

    if (opts->kind & MATLAB) {
      concat4(outfname,"A",opts->suffix,".mat");
      /* SRW -- this is binary data */
      if ( (fp = fopen(outfname,"wb")) == NULL) {
	printf("Couldn't open file\n");
	exit(1);
      }
      dump_M_to_matlab(fp, &Alist[1], num_real_nodes - 1, nnz - nnz0,"A");
      fclose(fp);
    }
  }

  /* save text files */
  if (opts->kind & TEXT && opts->dumpMats & DUMP_RL) {
    /* SRW -- this is ascii data */
    /*
    if ( (fp2 = fopen("M.dat","w")) == NULL) {
      printf("Couldn't open file\n");
      exit(1);
    }

    for(i = 0; i < num_fils; i++)
      Mrow[i] = 0;

    for(i = 0; i < num_mesh; i++) {
      fillMrow(indsys->Mlist, i, Mrow);
      dumpVec_totextfile(fp2, Mrow, num_fils);
    }

    fclose(fp2);
    */

    concat4(outfname,"R",opts->suffix,".dat");
    /* SRW -- this is ascii data */
    if ( (fp2 = fopen(outfname,"w")) == NULL) {
      printf("Couldn't open file\n");
      exit(1);
    }

    dumpVec_totextfile(fp2, R, num_fils);

    fclose(fp2);

    if (!indsys->dont_form_Z && indsys->opts->mat_vect_prod == DIRECT) {
      /* do imaginary part of Z */

      concat4(outfname,"L",opts->suffix,".dat");
      /* SRW -- this is ascii data */
      if ( (fp2 = fopen(outfname,"w")) == NULL) {
	printf("Couldn't open file\n");
	exit(1);
      }

      dumpMat_totextfile(fp2, Z, num_fils, num_fils);

      fclose(fp2);

    }
    else {
      if (indsys->dont_form_Z)
        printf("L matrix not dumped since L = 0 since fmin = 0\n");
      else
        printf("L matrix not dumped.  Run with mat_vect_prod = DIRECT if this is desired\n");
    }
  }

  /* Save Matlab files */

  if (opts->kind & MATLAB && opts->dumpMats & DUMP_RL) {
    concat4(outfname,"RL",opts->suffix,".mat");
    /* SRW -- this is binary data */
    if ( (fp = fopen(outfname,"wb")) == NULL) {
      printf("Couldn't open file\n");
      exit(1);
    }

    machine = 1000;
#ifdef DEC
    machine = 0000;
#endif

    dumbx =  (double *)malloc(num_fils*sizeof(double));
    dumby =  (double *)malloc(num_fils*sizeof(double));
    dumbz =  (double *)malloc(num_fils*sizeof(double));
    dumb =  (double *)malloc(num_fils*sizeof(double));
    if (dumb == NULL) {
      printf("no space for R\n");
      exit(1);
    }

    /* save and print matrices */

    /*
    for(i = 0; i < num_fils; i++)
      Mrow[i] = 0;

    for(i = 0; i < num_mesh; i++) {
      fillMrow(indsys->Mlist, i, Mrow);
      savemat_mod(fp, machine+100, "M", num_mesh, num_fils, 0, Mrow,
		  (double *)NULL, i, num_fils);
    }
    */

    /* this only saves the real part (savemat_mod.c modified) */
    savemat_mod(fp, machine+100, "R", 1, num_fils, 0, R, (double *)NULL,
		0, num_fils);

    if (!indsys->dont_form_Z && indsys->opts->mat_vect_prod == DIRECT) {
      /* do imaginary part of Z */
      for(i = 0; i < num_fils; i++) {
	savemat_mod(fp, machine+100, "L", num_fils, num_fils, 0, Z[i],
		    (double *)NULL, i, num_fils);
      }
    }
    else {
      if (indsys->dont_form_Z)
        printf("L matrix not dumped since L = 0 since fmin = 0\n");
      else
        printf("L matrix not dumped.  Run with mat_vect_prod = DIRECT if this is desired\n");
    }

    /* save vector of filament areas */
    for(tmps = indsys->segment; tmps != NULL; tmps = tmps->next) {
      for(j = 0; j < tmps->num_fils; j++) {
	tmpf = &(tmps->filaments[j]);
	dumb[tmpf->filnumber] = tmpf->area;
	dumbx[tmpf->filnumber] = tmpf->x[0];
	dumby[tmpf->filnumber] = tmpf->y[0];
	dumbz[tmpf->filnumber] = tmpf->z[0];
      }
    }
    savemat_mod(fp, machine+0, "areas", num_fils, 1, 0, dumb, (double *)NULL,0,
		num_fils);
    savemat_mod(fp, machine+0, "pos", num_fils, 3, 0, dumbx, (double *)NULL, 0,
		num_fils);
    savemat_mod(fp, machine+0, "pos", num_fils, 3, 0, dumby, (double *)NULL, 1,
		num_fils);
    savemat_mod(fp, machine+0, "pos", num_fils, 3, 0, dumbz, (double *)NULL, 1,
		num_fils);

    free(dumb);
    free(dumbx);
    free(dumby);
    free(dumbz);

    fclose(fp);
  }
}

void savecmplx(FILE *fp, char *name, CX **Z, int rows, int cols)
{
  int i,j;
  int machine;

  if (fp==NULL)
  {
    return;
  }

  machine = 1000;
#ifdef DEC
  machine = 0000;
#endif

  /* this only saves the real part (savemat_mod.c modified) one byte per call*/
  for(i = 0; i < rows; i++)
    for(j = 0; j < cols; j++)
      savemat_mod(fp, machine+100, name, rows, cols, 1, &Z[i][j].real,
		  (double *)NULL, i+j, 1);

  /* do imaginary part of Z */
  for(i = 0; i < rows; i++)
    for(j = 0; j < cols; j++)
      savemat_mod(fp, machine+100, name, rows, cols, 0, &Z[i][j].imag,
		  (double *)NULL, 1, 1);
}

/* saves a complex matrix more efficiently? */
void savecmplx2(FILE *fp, char *name, CX **Z, int rows, int cols)
{
  int i,j;
  int machine;
  static double *temp;
  static int colmax = 0;

  if (colmax < cols) {
    temp = (double *)malloc(cols*sizeof(double));
    colmax = cols;
  }

  machine = 1000;
#ifdef DEC
  machine = 0000;
#endif

  /* this only saves the real part (savemat_mod.c modified) one byte per call*/
  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++)
      temp[j] = Z[i][j].real;
    savemat_mod(fp, machine+100, name, rows, cols, 1, temp,
		(double *)NULL, i, cols);
  }

  /* do imaginary part of Z */
  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++)
      temp[j] = Z[i][j].imag;
    savemat_mod(fp, machine+100, name, rows, cols, 0, temp,
		(double *)NULL, 1, cols);
  }
}

/* This computes the product M*Z*Mt in a better way than oldformMZMt */
void formMZMt(SYS *indsys)
{
  int m,n,p;
  double tempsum, tempR, tempsumR;
  static double *tcol = NULL;   /* temporary storage for extra rows */
  int i,j, k, mesh, mesh2, nodeindx;
  int nfils, nmesh;
  MELEMENT *melem, *melem2;
  MELEMENT *mt, *mt2;    /* elements of M transpose */
  double **M, **L, *R;
  CX **Zm, *tempZ;
  int rows, cols, num_mesh;
  MELEMENT **Mlist, **Mt;

  Zm = indsys->MtZM;
  L = indsys->Z;
  R = indsys->R;
  M = indsys->M;
  nfils = rows = indsys->num_fils;
  nmesh = cols = num_mesh = indsys->num_mesh;
  Mlist = indsys->Mlist;
  Mt = indsys->Mtrans;

  if (nmesh > nfils) {
    fprintf(stderr, "Uh oh, more meshes than filaments, I'm confused\n");
    exit(1);
  }

  /* allocate a temporary column for Z*Mt */
  if (tcol == NULL)
    tcol = (double *)MattAlloc(rows, sizeof(double));

  /* this does L*(Mt)_i where (Mt)_i is a single column (row) of Mt (M) and
     saves it in the temp space, tcol */
  for(mesh = 0; mesh < num_mesh; mesh++) {
    for(j = 0; j< nfils; j++)
      tcol[j] = 0;
    /* note, this next set of nested loops could be reversed, but I think
       this might be more efficient for pipelining, ? */
    for(melem = Mlist[mesh]; melem != NULL; melem = melem->mnext)
      for(j = 0; j < nfils; j++)
	tcol[j] += L[j][melem->filindex]*melem->sign;
    for(mesh2 = 0; mesh2 < num_mesh; mesh2++) {
      tempsum = 0;
      for(melem2 = Mlist[mesh2]; melem2 != NULL; melem2 = melem2->mnext)
	tempsum += melem2->sign*tcol[melem2->filindex];
      Zm[mesh2][mesh].imag = tempsum;
    }
  }

  /* R is diagonal, so M*R*Mt is easy */
  for(i = 0; i < num_mesh; i++)
    for(j = 0; j < num_mesh; j++)
      Zm[i][j].real = 0;

  for(i = 0; i < nfils; i++)
    for(mt = Mt[i]; mt != NULL; mt = mt->mnext) {
      tempR = R[i]*mt->sign;
      for(mt2 = Mt[i]; mt2 != NULL; mt2 = mt2->mnext)
	Zm[mt2->filindex][mt->filindex].real += mt2->sign*tempR;
    }
}

void oldformMZMt(SYS *indsys)
{
  int m,n,p;
  double tempsum;
  static CX **trows = NULL;   /* temporary storage for extra rows */
  int i,j, k, mesh, nodeindx;
  int nfils, nmesh;
  MELEMENT *melem, *melem2;
  double **M, **L, *R;
  CX **Zm;
  int rows, cols, num_mesh;
  MELEMENT **Mlist;

  Zm = indsys->MtZM;
  L = indsys->Z;
  R = indsys->R;
  M = indsys->M;
  nfils = rows = indsys->num_fils;
  nmesh = cols = num_mesh = indsys->num_mesh;
  Mlist = indsys->Mlist;

  if (nmesh > nfils) {
    fprintf(stderr, "Uh oh, more meshes than filaments, I'm confused\n");
    exit(1);
  }

  /* allocate some extra rows so we can use Zm as temp space */
  if (rows > cols && trows == NULL)
    trows = (CX **)MatrixAlloc(rows - cols, cols, sizeof(CX));

  /* this does L*Mt and saves it in Zm.real.  This could be done better i
     think (use the fact that MZMt is symmetric)
     Also, could only track through each mesh list once, and store
     temporary values in Zm.real as we go along. (someday) */
  for(mesh = 0; mesh < num_mesh; mesh++) {
    for(j = 0; j < nfils; j++) {
      tempsum = 0;
      for(melem = Mlist[mesh]; melem != NULL; melem = melem->mnext)
	tempsum += L[j][melem->filindex]*melem->sign;
      if (j < nmesh)
	Zm[j][mesh].real = tempsum;
      else
	trows[j - nmesh][mesh].real = tempsum;
    }
  }

/*      savecmplx(fp, "step1", Zm, num_mesh, num_mesh); */

  /* this does M*(Zm.real) where Zm.real is hopefully L*Mt */
  for(mesh = 0; mesh < num_mesh; mesh++) {
    for(j = 0; j < nmesh; j++) {
      tempsum = 0;
      for(melem = Mlist[mesh]; melem != NULL; melem = melem->mnext) {
	if (melem->filindex < nmesh)
	  tempsum += Zm[melem->filindex][j].real*melem->sign;
	else
	  tempsum += trows[melem->filindex - nmesh][j].real*melem->sign;
      }
      Zm[mesh][j].imag = tempsum;
    }
  }
/*      savecmplx(fp, "step2", Zm, num_mesh, num_mesh); */

  /* R is diagonal, so M*R*Mt is easy */
  for(i = 0; i < num_mesh; i++) {
    for(j = i; j < num_mesh; j++) {
      tempsum = 0;
      /* brute force search for elements */
      for(melem = Mlist[i]; melem != NULL; melem = melem->mnext) {
	for(melem2 = Mlist[j]; melem2 != NULL; melem2 = melem2->mnext) {
	  if (melem2->filindex == melem->filindex)
	    tempsum += melem->sign*R[melem->filindex]*melem2->sign;
	}
      }
      Zm[i][j].real = Zm[j][i].real = tempsum;
    }
  }
  /* don't free the temp space */
  /*
  if (rows > cols) {
    for(i = 0; i < rows - cols; i++)
      free(trows[i]);
    free(trows);
  }
  */
}


/* forms the transpose of M.  Makes a linked list of each row.  Mlist is
    a linked list of its rows. */
/* Note: This uses the same struct as Mlist but in reality, each linked list
   is composed of mesh indices, not fil indices. (filindex is a mesh index) */
void formMtrans(SYS *indsys)
{
  int i, j, count;
  MELEMENT *m, *mt, *mt2, mtdum;
  MELEMENT **Mlist, **Mtrans, **Mtrans_check;
  int meshes, fils;
  int last_filindex;

  fils = indsys->num_fils;
  meshes = indsys->num_mesh;
  Mlist = indsys->Mlist;
  Mtrans = indsys->Mtrans;

  /* clear it (should be garbage only) */
  for(i = 0; i < fils; i++)
    Mtrans[i] = NULL;

  for(j = 0; j < meshes; j++)
  {
    last_filindex = -1;
    for(m = Mlist[j]; m != NULL; m = m->mnext)
    {
      mt = make_melement(indsys,j, NULL, m->sign);
      Mtrans[m->filindex] = insert_in_list(indsys,mt,Mtrans[m->filindex]);
      if (last_filindex == m->filindex)
	      printf("Internal Warning: Mesh %d contains the same filament multiple times\n",j);
      last_filindex = m->filindex;
    }
  }

#if 1==0
  /* this was the old inefficient way to from Mt */

  /* allocated for comparison */
  Mtrans_check = (MELEMENT **)MattAlloc(fils,sizeof(MELEMENT *));

  for(i = 0; i < fils; i++) {
    mt = &mtdum;
    /* scan through mesh list j looking for a filament i. (j,i) in M */
    for(j = 0; j < meshes; j++) {
      for(m = Mlist[j]; m->filindex < i && m->mnext != NULL; m = m->mnext)
	;
      if (m->filindex == i) {
	mt->mnext = make_melement(indsys,j, NULL, m->sign);
	mt = mt->mnext;
	if (m->mnext != NULL && m->mnext->filindex == i) {
	  for(count = 1; m->mnext->filindex == i; count++) {
	    m = m->mnext;
	    mt->mnext = make_melement(indsys,j, NULL, m->sign);
	    mt = mt->mnext;
	  }
	  printf("Internal Warning: Mesh %d contains the same filament %d times\n",j,count);
	}
      }
    }
    mt->mnext = NULL;
    Mtrans_check[i] = mtdum.mnext;
  }

  /* check old and new ways */
  for(i = 0; i < fils; i++) {
    for(mt = Mtrans[i], mt2 = Mtrans_check[i]; mt != NULL && mt2 != NULL;
                  mt = mt->mnext, mt2 = mt2->mnext)
      if (mt->filindex != mt2->filindex || mt->sign != mt2->sign)
        printf("different: %d  %d %d\n",i,mt->filindex, mt2->filindex);
    if (mt != mt2)
      printf("both not NULL:  mt: %u  mt2: %u\n", mt, mt2);
  }
#endif

}

void compare_meshes(MELEMENT *mesh1, MELEMENT *mesh2)
{

  while(mesh1 != NULL && mesh2 != NULL && mesh1->filindex == mesh2->filindex && mesh1->sign == mesh2->sign) {
    mesh1 = mesh1->mnext;
    mesh2 = mesh2->mnext;
  }

  if (mesh1 != NULL || mesh2 != NULL) {
    fprintf(stderr, "meshes don't match!\n");
    exit(1);
  }
}

void cx_dumpMat_totextfile(FILE *fp, CX **Z, int rows, int cols)
{
  int i, j;

  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++)
      fprintf(fp, "%13.6lg %+13.6lgj ", Z[i][j].real, Z[i][j].imag);
    fprintf(fp, "\n");
  }
  return;
}

void dumpMat_totextfile(FILE *fp, double **A, int rows, int cols)
{
  int i, j;

  for(i = 0; i < rows; i++) {
    for(j = 0; j < cols; j++)
      fprintf(fp, "%13.6lg ", A[i][j]);
    fprintf(fp, "\n");
  }
  return;
}

void dumpVec_totextfile(FILE *fp2, double *Vec, int size)
{
  dumpMat_totextfile(fp2, &Vec, 1, size);
}

void fillMrow(MELEMENT **Mlist, int mesh, double *Mrow)
{
  int i;
  MELEMENT *melem;

  if (mesh != 0)
    /* remove non-zeros from last call */
    for(melem = Mlist[mesh - 1]; melem != NULL; melem=melem->mnext)
      Mrow[melem->filindex] = 0;

  for(melem = Mlist[mesh]; melem != NULL; melem = melem->mnext)
    Mrow[melem->filindex] = melem->sign;
}

void dump_to_Ycond(FILE *fp, int cond, SYS *indsys)
{
  static char fname[20], tempstr[5];
  int maxiters = indsys->opts->maxiters;

  sprintf(tempstr, "%d", cond);

  strcpy(fname,"Ycond");
  strcat(fname,tempstr);
  savecmplx(fp, fname, indsys->FinalY, indsys->num_extern,
	    indsys->num_sub_extern);

  strcpy(fname,"resids");
  strcat(fname,tempstr);
  saveCarray(fp, fname, indsys->resids, indsys->num_sub_extern, maxiters);

  strcpy(fname,"resid_real");
  strcat(fname,tempstr);
  saveCarray(fp, fname, indsys->resid_real, indsys->num_sub_extern, maxiters);

  strcpy(fname,"resid_imag");
  strcat(fname,tempstr);
  saveCarray(fp, fname, indsys->resid_imag, indsys->num_sub_extern, maxiters);

  strcpy(fname,"niters");
  strcat(fname,tempstr);
  saveCarray(fp, fname, &(indsys->niters), 1, indsys->num_sub_extern);

}

void saveCarray(FILE *fp, char *fname, double **Arr, int rows, int cols)
{
  int i;
  int machine;

    machine = 1000;
#ifdef DEC
    machine = 0000;
#endif

  for(i = 0; i < rows; i++) {
    savemat_mod(fp, machine+100, fname, rows, cols, 0, Arr[i], (double *)NULL,
		i, cols);
  }
}

int nnz_inM(MELEMENT **Mlist, int num_mesh)
{
  int counter, i;
  MELEMENT *mesh;

  counter = 0;

  for(i = 0; i < num_mesh; i++)
    for(mesh = Mlist[i]; mesh != NULL; mesh = mesh->mnext)
      counter++;

  return counter;
}

void dump_M_to_text(FILE *fp, MELEMENT **Mlist, int num_mesh, int nnz)
{
  int counter, i;
  MELEMENT *mesh;

  counter = 0;

  for(i = 0; i < num_mesh; i++)
    for(mesh = Mlist[i]; mesh != NULL; mesh = mesh->mnext) {
      fprintf(fp, "%d\t%d\t%d\n", i+1, mesh->filindex+1, mesh->sign);
      counter++;
    }

  if (counter != nnz)
    fprintf(stderr,"Internal Warning: nnz %d != counter %d\n",nnz,counter);

}

void dump_M_to_matlab(FILE *fp, MELEMENT **Mlist, int num_mesh, int nnz,
    char *mname)
{

  double onerow[3];
  int i, counter;
  MELEMENT *mesh;

  counter = 0;

  for(i = 0; i < num_mesh; i++) {
    onerow[0] = i + 1;
    for(mesh = Mlist[i]; mesh != NULL; mesh = mesh->mnext) {
      onerow[1] = mesh->filindex + 1;
      onerow[2] = mesh->sign;
      savemat_mod(fp, machine+100, mname, nnz, 3, 0, onerow,
		  (double *)NULL, counter++, 3);
    }
  }

  if (counter != nnz)
    fprintf(stderr,"Internal Warning: nnz %d != counter %d\n",nnz,counter);
}

/* this picks one node in each tree to be a ground node */
void pick_ground_nodes(SYS *indsys)
{
  TREE *atree;
  SEGMENT *seg;
  seg_ptr segp;
  char type;
  NODES *node;

  for(atree = indsys->trees; atree != NULL; atree = atree->next) {
    type = atree->loops->path->seg.type;
    if (type == NORMAL)
      node = getrealnode(((SEGMENT *)atree->loops->path->seg.segp)->node[0]);
    else
      node = getrealnode(((PSEUDO_SEG *)atree->loops->path->seg.segp)->node[0]);

    if (node->index != -1) {
      fprintf(stderr,"huh? index == %d in pick_ground_node\n",node->index);
      exit(1);
    }
    node->index = 0;
  }
}

int pick_subset(strlist *portlist, SYS *indsys)
{
  strlist *oneport;
  EXTERNAL *ext;
  int counter;

  if (portlist == NULL) {
    for(ext = indsys->externals; ext != NULL; ext=ext->next)
      ext->col_Yindex = ext->Yindex;
    return indsys->num_extern;
  }

  for(ext = indsys->externals; ext != NULL; ext=ext->next)
    ext->col_Yindex = -1;

  counter = 0;
  for(oneport = portlist; oneport != NULL; oneport = oneport->next) {
    ext = get_external_from_portname(oneport->str,indsys);
    if (ext == NULL) {
      fprintf(stderr,"Error: unknown portname %s\n",oneport->str);
      exit(1);
    }

    ext->col_Yindex = counter;
    counter++;
  }

  return counter;
}

/* concatenates so that s1 = s1 + s2 + s3 + s4 */
void concat4(char *s1, char *s2, char *s3, char *s4)
{
  s1[0] = '\0';
  strcat(s1,s2);
  strcat(s1,s3);
  strcat(s1,s4);
}


double **sysMatrixAlloc(SYS* indsys, int rows, int cols, int size)
{
  double **temp;
  int i;

  temp=NULL;
  sysALLOC(temp,rows,double*,ON,IND,indsys,sysAllocTypePtrArray);
  if (temp == NULL)
  {
    printf("not enough space for matrix allocation\n");
    exit(1);
  }

  for(i = 0; i < rows; i++)
  {
    sysALLOC(temp[i],cols*size/sizeof(char),char,ON,IND,indsys,sysAllocTypeGeneric);
  }

  if (temp[rows - 1] == NULL)
  {
    printf("not enough space for matrix allocation\n");
    exit(1);
  }
  return(temp);
}


void sysDestroyMatrix (SYS* indsys, double*** ptr, int rows)
{
    int i;

    /* No pointer provided , so nothing to do */
    if (ptr==NULL)
    {
      return;
    }

   /* Link list isn't there, nothing to do */
   if (*ptr==NULL)
   {
     return;
   }

   /* Copy content */
   for(i = 0; i < rows; i++)
   {
     sysFree(indsys, ((*ptr)[i]));
   }

   sysFree(indsys, *ptr);
}

/** **/
/** Artifical memory management includes hash table and bulk allocation **/
/** by Oliver Wackerl **/
/** **/

/*
 Allocates bulk of specific data type if required and return one entry
 Uses a list of empty elements
 The chain pointer is required to be located in a  pointer field of regular content.
 Otherwise transfer with MPI transfer function breaks empty chain!

 There are multiple data types accelerated with bulk allocation and empty element chains.
 Wrote a single code for it, but it got a bit hard to read...
*/
void sysAllocateRecord ( void** PNTR, SYS* sys, void** EmptyList, unsigned long cnt, unsigned long typesize, unsigned int AllocatedType, void* baseaddr, void* nextaddr )
{
    ASSERT (PNTR!=NULL)
    ASSERT (sys!=NULL)
    ASSERT (EmptyList!=NULL)

    if (*EmptyList==NULL)
    {
      unsigned int  tcnt;
      unsigned long tsize;
      char*  ListPtr, *ListPtr1;

      tsize=cnt*typesize;

      /* Allocate block of records for allocation list. */
      ListPtr=NULL;
      CALLOC(ListPtr, tsize/sizeof(char), char, ON, AMSC);
      if (ListPtr == NULL)
      {
        fprintf(stderr, "Out of memory\n");
        return;
      }

      ListPtr1=ListPtr;
      /* Link all other new entries in Empty entry chain */
      for (tcnt=0; tcnt < cnt; tcnt++)
      {
        /* some pointer acrobatics to fill next pointer field of data type */
        unsigned long offset;
        char* tempc;
        void** tempv;
        offset=((char*)nextaddr)-((char*)baseaddr);

        tempc=(char*)ListPtr;
        tempc=&tempc[offset];
        tempv=(void**)tempc;

        *tempv= *EmptyList;
        *EmptyList = (void*)ListPtr;
        ListPtr=&ListPtr[typesize/sizeof(char)];
      }

      /* Register allocation */
      sysRecordAllocationFunc(sys, ListPtr1, tsize, AllocatedType, __FILE__, __LINE__ );
    }

    /* Now there must be empty elements available */
    ASSERT (*EmptyList!=NULL)

    /* Take one element of the empty chain */
    {
      /* some pointer acrobatics to address next pointer field of varying data type */
      unsigned long offset;
      char* tempc;
      void** tempv;
      offset=((char*)nextaddr)-((char*)baseaddr);

      *PNTR=*EmptyList;

      tempc=(char*)(*PNTR);
      tempc=&tempc[offset];
      tempv=(void**)tempc;

      *EmptyList= *tempv;
      /* Erase old pointer to next empty element */
      *tempv=NULL;
    }

}

void sysFreeRecord (void** PNTR, unsigned long typesize, void** EmptyList, void* baseaddr, void* nextaddr)
{
        /* some pointer acrobatics to fill next pointer field of varying data type */
        unsigned long offset;
        char* tempc;
        void** tempv;
        offset=((char*)nextaddr)-((char*)baseaddr);

         memset (*PNTR, 0, typesize);

        tempc=(char*)*PNTR;
        tempc=&tempc[offset];
        tempv=(void**)tempc;

        *tempv= *EmptyList;
        *EmptyList = *PNTR;
}

void sysAllocateBlockOfAllocationHashList( SYS* sys )
{
    unsigned int  tcnt;
    sysAllocationListHashPtr  ListPtr, ListPtr1;

    /* Begin `mulAllocateBlockOfAllocationList'. */
    ASSERT (sys!=NULL)

    /* Allocate block of records for allocation list. */
    ListPtr=NULL;
    CALLOC(ListPtr, (sysELEMENTS_PER_ALLOCATION+1), struct sysAllocationRecordHash, ON, AMSC);
    if (ListPtr == NULL)
    {

        return;
    }

    ListPtr1=ListPtr;
    /* Link all other new entries in EmptyAllocHashRecord chain */
    for (tcnt=0; tcnt < sysELEMENTS_PER_ALLOCATION+1; tcnt++)
    {
         ListPtr->NextRecord = sys->EmptyAllocHashRecord;
         sys->EmptyAllocHashRecord= ListPtr;
         ListPtr++;
    }

    /* Register istself */
    sysRecordAllocationFunc(sys, (char*) ListPtr1, (sysELEMENTS_PER_ALLOCATION+1)*sizeof(struct sysAllocationRecordHash), sysAllocTypeAllocationHashRecord, __FILE__, __LINE__ );
}

void sysAllocateBlockOfAllocationList( SYS* sys )
{
    unsigned int  tcnt;
    sysAllocationListPtr  ListPtr, ListPtr1;

    /* Begin `mulAllocateBlockOfAllocationList'. */
    ASSERT (sys!=NULL)

    /* Allocate block of records for allocation list. */
    ListPtr=NULL;
    CALLOC(ListPtr, (sysELEMENTS_PER_ALLOCATION+1), struct sysAllocationRecord, ON, AMSC);
    if (ListPtr == NULL)
    {

        return;
    }

    ListPtr1=ListPtr;
    /* Link all other new entries in EmptyAllocHashRecord chain */
    for (tcnt=0; tcnt < sysELEMENTS_PER_ALLOCATION+1; tcnt++)
    {
         ListPtr->NextRecord = sys->EmptyAllocRecord;
         sys->EmptyAllocRecord= ListPtr;
         ListPtr++;
    }

    /* Register itself */
    sysRecordAllocationFunc(sys, (char*) ListPtr1, (sysELEMENTS_PER_ALLOCATION+1)*sizeof(struct sysAllocationRecord), sysAllocTypeAllocationRecord, __FILE__, __LINE__ );
}


void sysRecordAllocationFunc( SYS* sys, char *AllocatedPtr, unsigned long AllocatedSize, unsigned int AllocatedType, const char* filename, unsigned int linenumber )
{
    #define bittabsize MAX((1u<<sysCalcHashWidth)/8,1)
    unsigned int  i, tsize, tsize1, done;
    unsigned char bittab[bittabsize];
    unsigned int thash;
    sysAllocationListHashPtr tblockLast;
    sysAllocationListPtr tallocblock;

    /*
     * If Allocated pointer is NULL, assume that malloc returned a NULL pointer,
    */
    ASSERT (sys!=NULL)

    ASSERT (AllocatedPtr!=NULL)
    ASSERT (AllocatedSize!=0)

    if ((AllocatedPtr == NULL) || (AllocatedSize==0))
    {
        return;
    }

    /* Replenish empty allocation record blocks if all are gone */
    if (sys->EmptyAllocRecord==NULL)
    {
       sysAllocateBlockOfAllocationList( sys );
    }
    ASSERT (sys->EmptyAllocRecord!=NULL)

    /* Grab one block from Empty Alloc Record list */
    tallocblock=sys->EmptyAllocRecord;
    sys->EmptyAllocRecord=tallocblock->NextRecord;

    /* Chain in in Allocation List */
    if (sys->TopOfAllocationList!=NULL)
    {
       sys->TopOfAllocationList->PrevRecord=tallocblock;
    }
    tallocblock->NextRecord=sys->TopOfAllocationList;
    tallocblock->PrevRecord=NULL;
    sys->TopOfAllocationList=tallocblock;

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
      thash=sysCalcHash( &AllocatedPtr[tsize1/sizeof(char)]  );

      /* Calculate bit position in bit map */
      bitaddr=thash>>3;
      bitmask=1<<(thash&7);

      /* If hash didn't occur so far add it */
      if (( bittab[ bitaddr] & bitmask) ==0)
      {
       sysAllocationListHashPtr tblock;

       /* Replenish empty hash record blocks if all are gone */
       if (sys->EmptyAllocHashRecord==NULL)
       {
         sysAllocateBlockOfAllocationHashList( sys );
       }
       ASSERT (sys->EmptyAllocHashRecord!=NULL)

       /* Remove from EmptyAllocHashRecord chain */
       tblock=sys->EmptyAllocHashRecord;
       sys->EmptyAllocHashRecord=tblock->NextRecord;

       /* Insert record in double linked list of specific hash chain */
       tblock->NextRecord = sys->TopOfAllocationHashList[thash];
       tblock->PrevRecord = NULL;
       if (sys->TopOfAllocationHashList[thash]!=NULL)
       {
          sys->TopOfAllocationHashList[thash]->PrevRecord=tblock;
       }
       sys->TopOfAllocationHashList[thash] = tblock;
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
        tsize1+=(1<<sysCalcHashGridSize);
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

void sysRemoveRecordAllocationFunc (SYS* sys, char *AllocatedPtr, const char* filename, unsigned int linenumber)
{
   /* Begin `mulRemoveRecordAllocationFunc'. */
   ASSERT (sys!=NULL)
   sysAllocationListHashPtr ListPtr;
   sysAllocationListPtr tallocblock;

   if (AllocatedPtr == NULL)
    {
        return;
    }

    ListPtr=sysFindRecordAllocation (sys, AllocatedPtr );
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
      sysAllocationListHashPtr EListPtr;

      /* Just to be sure no memory curruption is ongoing */
      ASSERT(ListPtr->AllocRecord==tallocblock)

      /* Remove element from double link chained table */
      if (ListPtr->PrevRecord==NULL)
      {
        sys->TopOfAllocationHashList[ListPtr->hashtable]=ListPtr->NextRecord;
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
      memset(EListPtr, 0, sizeof(struct sysAllocationRecordHash));

      /* Recycle entry */
      EListPtr->NextRecord=sys->EmptyAllocHashRecord;
      sys->EmptyAllocHashRecord=EListPtr;
    }

    /* Remove allocation entry from chain*/
    if (tallocblock->PrevRecord==NULL)
    {
      sys->TopOfAllocationList=tallocblock->NextRecord;
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
    memset(tallocblock, 0, sizeof(struct sysAllocationRecord));

    /* Link it to empty chain */
    tallocblock->NextRecord=sys->EmptyAllocRecord;
    sys->EmptyAllocRecord=tallocblock;
}




sysAllocationListHashPtr sysFindRecordAllocation (SYS* sys, char *Ptr )
{
/* Begin `mulFindRecordAllocation'. */

  sysAllocationListHashPtr  ListPtr;
  ASSERT (sys!=NULL)
  unsigned int thash;
  unsigned int cnt;

    if (Ptr==NULL)
    {
        return NULL;
    }

    thash=sysCalcHash(Ptr);
    ListPtr = sys->TopOfAllocationHashList[thash];
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

void sysALLOC_func(void** PNTR,unsigned long NUM,unsigned long TYPESIZE,int FLAG, int MTYP, unsigned int ATYP, SYS* sys, const char* filename, unsigned int linenumber)
{
  unsigned long tsize;
  tsize=(NUM*TYPESIZE);

  ASSERT (sys!=NULL)

  if ((*PNTR)!=NULL)
  {
    free(*PNTR);
    sysRemoveRecordAllocationFunc(sys, *PNTR, filename, linenumber);
  }

  *PNTR=NULL;

  if (tsize==0)
  {
    return;
  }

  if (ATYP==sysAllocTypeMElement)
  {
    MELEMENT temp;
    ASSERT(NUM==1)
    sysAllocateRecord ( PNTR, sys, (void**)&sys->EmptyMElement, sysMELEMENTS_PER_ALLOCATION, sizeof(MELEMENT), sysAllocTypeMElement, &temp, &temp.mnext );
    return;
  }

  if (ATYP==sysAllocTypeSegList)
  {
    SEGLIST temp;
    ASSERT(NUM==1)
    sysAllocateRecord ( PNTR, sys, (void**)&sys->EmptySEGLIST, sysSEGLIST_PER_ALLOCATION, sizeof(SEGLIST), sysAllocTypeSegList, &temp, &temp.next );
    return;
  }

  if (ATYP==sysAllocTypeSegment)
  {
    SEGMENT temp;
    ASSERT(NUM==1)
    sysAllocateRecord ( PNTR, sys, (void**)&sys->EmptySEGMENT, sysSEGMENT_PER_ALLOCATION, sizeof(SEGMENT), sysAllocTypeSegment, &temp, &temp.next );
    return;
  }

/*
  if (ATYP==sysAllocTypeFilament)
  {
    ASSERT(NUM==1)
    sysAllocateRecord ( PNTR, sys, (void**)&sys->EmptyFILAMENT, sysFILAMENT_PER_ALLOCATION, sizeof(FILAMENT), sysAllocTypeFilament );
    return;
  }

  if (ATYP==sysAllocTypeCharge)
  {
    ASSERT(NUM==1)
    sysAllocateRecord ( PNTR, sys, (void**)&sys->EmptyCHARGE, sysCHARGE_PER_ALLOCATION, sizeof(charge), sysAllocTypeCharge );
    return;
  }
*/

  if (ATYP==sysAllocTypeSurface)
  {
    surface temp;
    ASSERT(NUM==1)
    sysAllocateRecord ( PNTR, sys, (void**)&sys->EmptySURFACE, sysSURFACE_PER_ALLOCATION, sizeof(surface), sysAllocTypeSurface, &temp, &temp.next );
    return;
  }

  if (ATYP==sysAllocTypeNode)
  {
    NODES temp;
    ASSERT(NUM==1)
    sysAllocateRecord ( PNTR, sys, (void**)&sys->EmptyNODE, sysNODE_PER_ALLOCATION, sizeof(NODES), sysAllocTypeNode, &temp, &temp.next );
    return;
  }

  if (ATYP==sysAllocTypeSpath)
  {
    SPATH temp;
    ASSERT(NUM==1)
    sysAllocateRecord ( PNTR, sys, (void**)&sys->EmptySPATH, sysSPATH_PER_ALLOCATION, sizeof(SPATH), sysAllocTypeSpath, &temp, &temp.next );
    return;
  }

  ALLOC_func(0, PNTR, NUM, TYPESIZE, FLAG, MTYP, filename, linenumber);
  sysRecordAllocationFunc(sys,(char*)*PNTR,tsize,ATYP,filename,linenumber);
}

void sysReallocFunc (SYS* indsys, void** pPtr, unsigned long Size, const char* filename, unsigned int linenumber )
{
    void *oldPtr;
    sysAllocationListHashPtr  ListPtr;
    unsigned int typ;

   ASSERT (indsys!=NULL)
   ASSERT (pPtr!=NULL)
   ASSERT (*pPtr!=NULL)
   ASSERT (Size>0)

   oldPtr= *pPtr;

    ListPtr = sysFindRecordAllocation(indsys, oldPtr);
    ASSERT (ListPtr!=NULL)
    typ= ListPtr->AllocRecord->AllocatedType;

    (*pPtr) = realloc(*pPtr, Size);

    sysRemoveRecordAllocationFunc(indsys, (char*)oldPtr, filename, linenumber);

    sysRecordAllocationFunc(indsys, (char*)*pPtr, Size, typ, filename, linenumber);

}

void sysFreeFunc (SYS* indsys, void** pPtr, const char* filename, unsigned int linenumber )
{
   sysAllocationListHashPtr ListPtr;
   sysAllocationListHashPtr tblock;
   unsigned int  i, tsize, tsize1, done;
   unsigned char bittab[(1u<<sysCalcHashWidth)/8];
   char* ptr;

   ASSERT (indsys!=NULL)
   ASSERT (pPtr!=NULL)
   ASSERT (*pPtr!=NULL)

   ListPtr=sysFindRecordAllocation(indsys, *pPtr);
   ASSERT (ListPtr!=NULL)
   if (ListPtr->AllocRecord->AllocatedType==sysAllocTypeMElement)
   {
      MELEMENT temp;
      sysFreeRecord(pPtr, sizeof(MELEMENT), (void**)&indsys->EmptyMElement, &temp, &temp.mnext);
      return;
   }
   if (ListPtr->AllocRecord->AllocatedType==sysAllocTypeSegList)
   {
      SEGLIST temp;
      sysFreeRecord(pPtr, sizeof(SEGLIST), (void**)&indsys->EmptySEGLIST, &temp, &temp.next);
      return;
   }
   if (ListPtr->AllocRecord->AllocatedType==sysAllocTypeSurface)
   {
      surface temp;
      sysFreeRecord(pPtr, sizeof(surface), (void**)&indsys->EmptySURFACE, &temp, &temp.next);
      return;
   }
   if (ListPtr->AllocRecord->AllocatedType==sysAllocTypeSegment)
   {
     SEGMENT temp;
     sysFreeRecord(pPtr, sizeof(SEGMENT), (void**)&indsys->EmptySEGMENT, &temp, &temp.next);
     return;
   }
   if (ListPtr->AllocRecord->AllocatedType==sysAllocTypeNode)
   {
     NODES temp;
     sysFreeRecord(pPtr, sizeof(NODES), (void**)&indsys->EmptyNODE, &temp, &temp.next);
     return;
   }

   if (ListPtr->AllocRecord->AllocatedType==sysAllocTypeSpath)
   {
     SPATH temp;
     sysFreeRecord(pPtr, sizeof(SPATH), (void**)&indsys->EmptySPATH, &temp, &temp.next);
     return;
   }


   ptr=(char*)*pPtr;

   /* Free memory block */
   free(*pPtr);

   /* Delete link to block */
   (*pPtr) = NULL;

    sysRemoveRecordAllocationFunc (indsys, ptr,filename,linenumber);
}

unsigned long sysGetMemoryBlockLength (SYS* sys, char* AllocatedPtr)
{
   sysAllocationListHashPtr  ListPtr;
    ASSERT (sys!=NULL)

    if ((AllocatedPtr == NULL) || (sys==NULL))
    {
        return 0;
    }

    ListPtr = sysFindRecordAllocation(sys, AllocatedPtr);
    if (ListPtr!=NULL)
    {
        return ListPtr->AllocRecord->AllocatedSize;
    }

    return 0;
}

void sysDestroy( SYS** sys )

{
  sysAllocationListPtr  ListPtr;
  unsigned int i;

  ASSERT(sys!=NULL)
  if ((*sys)==NULL)
  {
    /* Nothing to do */
    return;
  }

  /* Make a copy of TopAllocationList entry */
  /* As soon block containing Matrix is gone */
  /* there is no access anymore */
  ListPtr=(*sys)->TopOfAllocationList;
  (*sys)=NULL;

  /* Sequentially step through the list of allocated pointers freeing pointers
   * along the way. */

  /* Run twice over list first run remove all data blocks */
  /* second run remove allocation records itself */
  for (i=0; i<2; i++)
  {
      unsigned int cnt;
      sysAllocationListPtr  tListPtr;

      cnt=0;
      tListPtr = ListPtr;
      while (tListPtr != NULL)
      {
          if ((i==1) || (tListPtr->AllocatedType!=sysAllocTypeAllocationRecord))
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
            memset(tptr1, 0, sizeof(struct sysAllocationRecord));

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



#define sysMarkAlloc(sys, ptr, flag) sysMarkAllocFunc(sys, (void*)ptr, flag, __FILE__, __LINE__)
int sysMarkAllocFunc (SYS* sys, void* ptr, unsigned int flag, const char* filename, unsigned int linenumber)
{
    ASSERT (sys!=NULL)

    sysAllocationListHashPtr  ListPtr;

    if (ptr==NULL)
      return 0;

    ListPtr = sysFindRecordAllocation(sys, (char*)ptr);
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

void sysAllocationAnalysis (SYS* sys)
{
    unsigned int allocs, tallocs, round;
    unsigned long allocsize;
    unsigned int found;
    sysAllocationListPtr  ListPtr;
    unsigned int population[100];

    if (sys==NULL)
    {
      return;
    }

    /* Reset Marks, mark all AllocationRecords as linked blocks */
    {
      unsigned int i;
      unsigned int tcnt;
      for (i=0; i<100; i++)
      {
        population[i]=0;
      }

      allocs=0;
      allocsize=0;
      ListPtr=sys->TopOfAllocationList;
      while (ListPtr!=NULL)
      {
        ASSERT(ListPtr->AllocatedSize>0)
        ASSERT(ListPtr->AllocatedPtr!=NULL)
       {
          allocs++;
          ListPtr->Mark=0;
          allocsize+=ListPtr->AllocatedSize;
          population[ListPtr->AllocatedType]++;
        }
        ListPtr=ListPtr->NextRecord;
      }

    }
    printf("Allocs %u, Size: %lu kB\n",allocs, allocsize>>10);
    {
      unsigned int i;
      for (i=0; i<100; i++)
      {
        if (population[i]>0)
        {
          printf("Typ %u: %u\n",i,population[i]);
        }
      }
    }

    /* mark block with root structure */
    sysMarkAlloc(sys,sys, 1);

    /* Run iterative over list until there isn't anything to mark anymore */
    found=1;
    tallocs=0;
    round=0;
    while (found>0)
    {
      unsigned int stat_Generic =0;
      unsigned int stat_PtrArray =0;
      unsigned int stat_AllocationRecord =0;
      unsigned int stat_indsys =0;
      unsigned int stat_indopts =0;
      unsigned int tcnt;

      printf("Round %u: %u / %u (%u%%) processed \r",round++,tallocs, allocs, (tallocs*100)/allocs); fflush(stdout);

      found=0;
      ListPtr=sys->TopOfAllocationList;
      while (ListPtr!=NULL)
      {

        if (ListPtr->Mark==1)
        {
            unsigned int err;
            found++;
            tallocs++;
            ListPtr->Mark=2;
            err=0;
            switch (ListPtr->AllocatedType)
            {
            case AllocTypeGeneric:
              stat_Generic++;
              /* Nothing to do */
            break;

            case sysAllocTypeAllocationHashRecord:
            {
             sysAllocationListHashPtr  pElement;
             unsigned int i,cnt;
             stat_AllocationRecord++;
             pElement=(sysAllocationListHashPtr)ListPtr->AllocatedPtr;
             cnt=ListPtr->AllocatedSize/sizeof(struct sysAllocationRecordHash);
               for (i=0; i<cnt; i++)
               {
                 unsigned int flag;
//                 err|=sysMarkAlloc(sys,pElement->AllocatedPtr,0);
                 err|=sysMarkAlloc(sys,pElement->NextRecord,1);
                 err|=sysMarkAlloc(sys,pElement->PrevRecord,1);
                 err|=sysMarkAlloc(sys,pElement->NextPeerRecord,1);
                 err|=sysMarkAlloc(sys,pElement->PrevPeerRecord,1);
                 err|=sysMarkAlloc(sys,pElement->AllocRecord,1);
                 pElement++;
               }

            }
            break;

            case sysAllocTypeAllocationRecord:
            {
             sysAllocationListPtr  pElement;
             unsigned int i,cnt;
             stat_AllocationRecord++;
             pElement=(sysAllocationListPtr)ListPtr->AllocatedPtr;
             cnt=ListPtr->AllocatedSize/sizeof(struct sysAllocationRecord);
               for (i=0; i<cnt; i++)
               {
                 err|=sysMarkAlloc(sys,pElement->AllocatedPtr,0);
                 err|=sysMarkAlloc(sys,pElement->NextRecord,1);
                 err|=sysMarkAlloc(sys,pElement->PrevRecord,1);
                 pElement++;
               }
            }
            break;

            case sysAllocTypePtrArray:
            {
              unsigned int tcount;
              int** tptr;
              unsigned int i;
              stat_PtrArray++;
              tcount=ListPtr->AllocatedSize/sizeof(int*);
              tptr=(int**)ListPtr->AllocatedPtr;
              for (i=0; i<tcount; i++)
              {
                err|=sysMarkAlloc(sys, (*tptr),1 );
                tptr++;
              }
            }
            break;

            case AllocTypessystem:
            { ssystem *tsys;
              tsys=(ssystem*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(ssystem));
              err|=sysMarkAlloc(sys,tsys->cond_names,1 );
              err|=sysMarkAlloc(sys,tsys->q,1);
              err|=sysMarkAlloc(sys,tsys->p,1);
              err|=sysMarkAlloc(sys,tsys->panels,1);
              err|=sysMarkAlloc(sys,tsys->cubes,1);
              err|=sysMarkAlloc(sys,tsys->multilist,1);
              err|=sysMarkAlloc(sys,tsys->locallist,1);
              err|=sysMarkAlloc(sys,tsys->directlist,1);
              err|=sysMarkAlloc(sys,tsys->precondlist,1);
              err|=sysMarkAlloc(sys,tsys->revprecondlist,1);
              err|=sysMarkAlloc(sys,tsys->is_dummy,1);
              err|=sysMarkAlloc(sys,tsys->is_dielec,1);
              err|=sysMarkAlloc(sys,tsys->indsys,1);
            }
            break;

            case AllocTypecube:
            {
              cube *tcube;
              tcube=(cube*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(cube));
              err|=sysMarkAlloc(sys,tcube->mnext,1);
              err|=sysMarkAlloc(sys,tcube->upnumeles,1);
              err|=sysMarkAlloc(sys,tcube->upvects,1);
              err|=sysMarkAlloc(sys,tcube->multi,1);
              err|=sysMarkAlloc(sys,tcube->upmats,1);
              err|=sysMarkAlloc(sys,tcube->is_dummy,1);
              err|=sysMarkAlloc(sys,tcube->is_dielec,1);
              err|=sysMarkAlloc(sys,tcube->lnext,1);
              err|=sysMarkAlloc(sys,tcube->downnumeles,1);
              err|=sysMarkAlloc(sys,tcube->downvects,1);
              err|=sysMarkAlloc(sys,tcube->local,1);
              err|=sysMarkAlloc(sys,tcube->downmats,1);
              err|=sysMarkAlloc(sys,tcube->interList,1);
              err|=sysMarkAlloc(sys,tcube->enext,1);
              err|=sysMarkAlloc(sys,tcube->evalnumeles,1);
              err|=sysMarkAlloc(sys,tcube->evalvects,1);
              err|=sysMarkAlloc(sys,tcube->eval,1);
              err|=sysMarkAlloc(sys,tcube->evalmats,1);
              err|=sysMarkAlloc(sys,tcube->eval_isQ2P,1);
              err|=sysMarkAlloc(sys,tcube->dnext,1);
              err|=sysMarkAlloc(sys,tcube->pnext,1);
              err|=sysMarkAlloc(sys,tcube->rpnext,1);
              err|=sysMarkAlloc(sys,tcube->directnumeles,1);
              err|=sysMarkAlloc(sys,tcube->directq,1);
              err|=sysMarkAlloc(sys,tcube->directmats,1);
              err|=sysMarkAlloc(sys,tcube->precondmats,1);
              err|=sysMarkAlloc(sys,tcube->directlu,1);
              err|=sysMarkAlloc(sys,tcube->precond,1);
              err|=sysMarkAlloc(sys,tcube->prevectq,1);
              err|=sysMarkAlloc(sys,tcube->prevectp,1);
              err|=sysMarkAlloc(sys,tcube->nbr_is_dummy,1);
              err|=sysMarkAlloc(sys,tcube->chgs,1);
              err|=sysMarkAlloc(sys,tcube->nbrs,1);
              err|=sysMarkAlloc(sys,tcube->kids,1);
              err|=sysMarkAlloc(sys,tcube->parent,1);
            }
            break;

            case AllocTypecharge:
            {
              charge *tcharge;
              tcharge=(charge*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(charge));
              err|=sysMarkAlloc(sys,tcharge->next,1);
              err|=sysMarkAlloc(sys,tcharge->multipole,1);
              err|=sysMarkAlloc(sys,tcharge->surf,1);
              err|=sysMarkAlloc(sys,tcharge->pos_dummy,1);
              err|=sysMarkAlloc(sys,tcharge->neg_dummy,1);
              err|=sysMarkAlloc(sys,tcharge->fil,1);
            }
            break;

            case AllocTypesurface:
            {
              surface *tsurface;
              tsurface=(surface*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(surface));
              err|=sysMarkAlloc(sys,tsurface->title,1);
              err|=sysMarkAlloc(sys,tsurface->name,1);
              err|=sysMarkAlloc(sys,tsurface->panels,1);
              err|=sysMarkAlloc(sys,tsurface->group_name,1);
              err|=sysMarkAlloc(sys,tsurface->next,1);
              err|=sysMarkAlloc(sys,tsurface->prev,1);
            }
            break;

            case AllocTypeName:
            {
              Name *tname;
              tname=(Name*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(Name));
              err|=sysMarkAlloc(sys,tname->name,1);
              err|=sysMarkAlloc(sys,tname->next,1);
              err|=sysMarkAlloc(sys,tname->alias_list,1);
            }
            break;


            case sysAllocTypeindsys:
            { SYS *tsys;
              unsigned int i;
              stat_indsys++;
              tsys=(SYS*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(SYS));
              err|=sysMarkAlloc(sys,tsys->planes,1);
              err|=sysMarkAlloc(sys,tsys->endplane,1);
              err|=sysMarkAlloc(sys,tsys->segment,1);
              err|=sysMarkAlloc(sys,tsys->endseg,1);
              err|=sysMarkAlloc(sys,tsys->nodes,1);
              err|=sysMarkAlloc(sys,tsys->endnode,1);
              err|=sysMarkAlloc(sys,tsys->pseudo_nodes,1);
              err|=sysMarkAlloc(sys,tsys->meshsect,1);
              err|=sysMarkAlloc(sys,tsys->diagL,1);
              /* err|=err|=sysMarkAlloc(sys,tsys->sparMatrix,1); */
              err|=sysMarkAlloc(sys,tsys->Alist,1);
              err|=sysMarkAlloc(sys,tsys->M,1);
              err|=sysMarkAlloc(sys,tsys->Mlist,1);
              err|=sysMarkAlloc(sys,tsys->m_info,1);
              err|=sysMarkAlloc(sys,tsys->Mtrans,1);
              err|=sysMarkAlloc(sys,tsys->Precond,1);
              err|=sysMarkAlloc(sys,tsys->trees,1);
              err|=sysMarkAlloc(sys,tsys->externals,1);
              err|=sysMarkAlloc(sys,tsys->Z,1);
              err|=sysMarkAlloc(sys,tsys->R,1);
              err|=sysMarkAlloc(sys,tsys->MtZM,1);
              err|=sysMarkAlloc(sys,tsys->FinalY,1);
              err|=sysMarkAlloc(sys,tsys->resids,1);
              err|=sysMarkAlloc(sys,tsys->resid_real,1);
              err|=sysMarkAlloc(sys,tsys->resid_imag,1);
              err|=sysMarkAlloc(sys,tsys->niters,1);
              err|=sysMarkAlloc(sys,tsys->opts,1);
              err|=sysMarkAlloc(sys,tsys->Ar,1);
              err|=sysMarkAlloc(sys,tsys->Br,1);
              err|=sysMarkAlloc(sys,tsys->Cr,1);
              err|=sysMarkAlloc(sys,tsys->Dr,1);
              err|=sysMarkAlloc(sys,tsys->title,1);
              err|=sysMarkAlloc(sys,tsys->chglist,1);
              err|=sysMarkAlloc(sys,tsys->sys,1);
              for (i=0; i<(1u<<sysCalcHashWidth); i++)
              {
                err|=sysMarkAlloc(sys,tsys->TopOfAllocationHashList[i],1);
              }
              err|=sysMarkAlloc(sys,tsys->EmptyAllocHashRecord,1);
              err|=sysMarkAlloc(sys,tsys->TopOfAllocationList,1);
              err|=sysMarkAlloc(sys,tsys->EmptyAllocRecord,1);
              err|=sysMarkAlloc(sys,tsys->EmptyMElement,1);
              err|=sysMarkAlloc(sys,tsys->EmptySEGLIST,1);
              err|=sysMarkAlloc(sys,tsys->EmptyFILAMENT,1);
              err|=sysMarkAlloc(sys,tsys->EmptyCHARGE,1);
              err|=sysMarkAlloc(sys,tsys->EmptySURFACE,1);
              err|=sysMarkAlloc(sys,tsys->EmptySEGMENT,1);
              err|=sysMarkAlloc(sys,tsys->EmptyNODE,1);
              err|=sysMarkAlloc(sys,tsys->EmptySPATH,1);
#ifdef SRW0814
              err|=sysMarkAlloc(sys,tsys->tab.entries,1);
#endif
            }
            break;

            case sysAllocTypeindopts:
            {
             ind_opts *tind_opts;
             stat_indopts++;
             tind_opts=(ind_opts*)ListPtr->AllocatedPtr;
             ASSERT (ListPtr->AllocatedSize==sizeof(ind_opts));
             err|=sysMarkAlloc(sys,tind_opts->portlist,1 );
             err|=sysMarkAlloc(sys,tind_opts->suffix,1 );
             err|=sysMarkAlloc(sys,tind_opts->fname,1 );
            }
            break;

            case sysAllocTypestrlist:
            {
             strlist *tstrlist;
             tstrlist=(strlist*)ListPtr->AllocatedPtr;
             ASSERT (ListPtr->AllocatedSize==sizeof(strlist));
             err|=sysMarkAlloc(sys,tstrlist->str,1 );
             err|=sysMarkAlloc(sys,tstrlist->next,1 );
            }
            break;

            case sysAllocTypeNode:
            {
              NODES* tnode;
              unsigned int i,tcnt;
              tnode=(NODES*)ListPtr->AllocatedPtr;
              tcnt=ListPtr->AllocatedSize/sizeof(NODES);
              for (i=0; i<tcnt; i++)
              {
                err|=sysMarkAlloc(sys,tnode->name,1 );
                err|=sysMarkAlloc(sys,tnode->equiv,1 );
                err|=sysMarkAlloc(sys,tnode->to_end,1 );
                err|=sysMarkAlloc(sys,tnode->connected_segs,1 );
                err|=sysMarkAlloc(sys,tnode->gp,1 );
                err|=sysMarkAlloc(sys,tnode->treeptr,1 );
                err|=sysMarkAlloc(sys,tnode->gp_node,1 );
                err|=sysMarkAlloc(sys,tnode->next,1 );
                tnode++;
              }
            }
            break;

            case sysAllocTypeSegment:
            {
              SEGMENT* tseg;
              unsigned tcount, tcnt;
              tseg=(SEGMENT*)ListPtr->AllocatedPtr;
              tcount=ListPtr->AllocatedSize/sizeof(SEGMENT);
              for (tcnt=0; tcnt<tcount; tcnt++)
              {
                err|=sysMarkAlloc(sys,tseg->name,1 );
                err|=sysMarkAlloc(sys,tseg->widthdir,1 );
                err|=sysMarkAlloc(sys,tseg->node[0],1 );
                err|=sysMarkAlloc(sys,tseg->node[1],1 );
                err|=sysMarkAlloc(sys,tseg->filaments,1 );
                err|=sysMarkAlloc(sys,tseg->loops,1 );
                err|=sysMarkAlloc(sys,tseg->next,1 );
              }
            }
            break;

            case sysAllocTypeGroundplane:
            {
              GROUNDPLANE *tgp;
              tgp=(GROUNDPLANE*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(GROUNDPLANE));
              err|=sysMarkAlloc(sys,tgp->name,1);
              err|=sysMarkAlloc(sys,tgp->grid1,1);
              err|=sysMarkAlloc(sys,tgp->grid2,1);
              err|=sysMarkAlloc(sys,tgp->segs1,1);
              err|=sysMarkAlloc(sys,tgp->segs2,1);
              err|=sysMarkAlloc(sys,tgp->pnodes,1);
              err|=sysMarkAlloc(sys,tgp->lastnode,1);
              err|=sysMarkAlloc(sys,tgp->innode,1);
              err|=sysMarkAlloc(sys,tgp->outnode,1);
              err|=sysMarkAlloc(sys,tgp->usernodes,1);
              err|=sysMarkAlloc(sys,tgp->fake_seg_list,1);
              err|=sysMarkAlloc(sys,tgp->list_of_holes,1);
              err|=sysMarkAlloc(sys,tgp->list_of_contacts,1);
              err|=sysMarkAlloc(sys,tgp->usernode_coords,1);
              err|=sysMarkAlloc(sys,tgp->filename,1);
              err|=sysMarkAlloc(sys,tgp->nonuni,1);
              err|=sysMarkAlloc(sys,tgp->indsys,1);
              err|=sysMarkAlloc(sys,tgp->next,1);
            }
            break;

            case sysAllocTypePseudoNode:
            {
              PSEUDO_NODE* tnode;
              tnode=(PSEUDO_NODE*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(PSEUDO_NODE));
              err|=sysMarkAlloc(sys,tnode->node,1);
              err|=sysMarkAlloc(sys,tnode->name,1);
              err|=sysMarkAlloc(sys,tnode->next,1);
            }
            break;

            case sysAllocTypeSpath:
            {
              SPATH* tspath;
              unsigned i, cnt;
              tspath=(SPATH*)ListPtr->AllocatedPtr;
              cnt=ListPtr->AllocatedSize/sizeof(SPATH);
              for (i=0; i<cnt; i++)
              {
               err|=sysMarkAlloc(sys,tspath->seg.segp,1);
               err|=sysMarkAlloc(sys,tspath->next,1);
              }
            }
            break;

            case sysAllocTypeMElement:
            {
              MELEMENT* tmelement;
              unsigned tcount, tcnt;
              tmelement=(MELEMENT*)ListPtr->AllocatedPtr;
              tcount=ListPtr->AllocatedSize/sizeof(MELEMENT);
              for (tcnt=0; tcnt<tcount; tcnt++)
              {
                err|=sysMarkAlloc(sys,tmelement->fil,1);
                err|=sysMarkAlloc(sys,tmelement->mnext,1);
                tmelement++;
              }
            }
            break;

            case sysAllocTypeContactList:
            {
              ContactList* tcontact;
              tcontact=(ContactList*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(ContactList));
              err|=sysMarkAlloc(sys,tcontact->func,1);
              err|=sysMarkAlloc(sys,tcontact->name,1);
              err|=sysMarkAlloc(sys,tcontact->vals,1);
              err|=sysMarkAlloc(sys,tcontact->next,1);
            }
            break;

            case sysAllocTypeBi:
            {
              Bi* tbi;
              tbi=(Bi*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(Bi));
              err|=sysMarkAlloc(sys,tbi->child1,1);
              err|=sysMarkAlloc(sys,tbi->child2,1);
            }
            break;

            case sysAllocTypeGrid_2d:
            {
              Grid_2d* tgrid_2d;
              tgrid_2d=(Grid_2d*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(Grid_2d));
              err|=sysMarkAlloc(sys,tgrid_2d->kids,1);
            }
            break;

            case sysAllocTypeInd_list:
            {
              int_list* tlist;
              tlist=(int_list*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(int_list));
              err|=sysMarkAlloc(sys,tlist->next,1);
            }
            break;

            case sysAllocTypeNonuni_gp:
            {
              Nonuni_gp *tgp;
              tgp=(Nonuni_gp*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(Nonuni_gp));
              err|=sysMarkAlloc(sys,tgp->root_cell,1);
              err|=sysMarkAlloc(sys,tgp->z_pts,1);
              err|=sysMarkAlloc(sys,tgp->thick,1);
              err|=sysMarkAlloc(sys,tgp->z_c,1);
              err|=sysMarkAlloc(sys,tgp->nodelist,1);
              err|=sysMarkAlloc(sys,tgp->grndp,1);
            }
            break;

            case sysAllocTypeG_edges:
            {
              G_edges *tedge;
              unsigned int i;
              tedge=(G_edges*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(G_edges));
              for (i=0; i<NUM_E_CELLS; i++)
              {
                err|=sysMarkAlloc(sys,tedge->cells[i],1);
              }
              err|=sysMarkAlloc(sys,tedge->children,1);
            }
            break;

            case sysAllocTypeGcell:
            {
              Gcell *tgcell;
              unsigned int i,cnt,j ;
              tgcell=(Gcell*)ListPtr->AllocatedPtr;
              cnt=ListPtr->AllocatedSize/sizeof(Gcell);
              for (j=0; j<cnt; j++)
              {
                err|=sysMarkAlloc(sys,tgcell->children,1);
                err|=sysMarkAlloc(sys,tgcell->parent,1);
                for (i=0; i<MAX(NUMEDGES,NUMNODES); i++)
                {
                  err|=sysMarkAlloc(sys,tgcell->bndry.edges[i],1);
                }
                tgcell++;
              }
            }
            break;

            case sysAllocTypeG_nodes:
            {
              G_nodes *tgnode;
              unsigned int i;
              tgnode=(G_nodes*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(G_nodes));
              for (i=0; i<NUM_N_CELLS; i++)
              {
                err|=sysMarkAlloc(sys,tgnode->cells[i],1);
              }
              for (i=0; i<NUMADJ; i++)
              {
                err|=sysMarkAlloc(sys,tgnode->adjacent[i],1);
              }
              err|=sysMarkAlloc(sys,tgnode->e_segs,1);
              err|=sysMarkAlloc(sys,tgnode->n_segs,1);
              err|=sysMarkAlloc(sys,tgnode->prev,1);
              err|=sysMarkAlloc(sys,tgnode->next,1);
            }
            break;

            case sysAllocTypeLlist:
            {
              Llist *tlist;
              tlist=(Llist*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(Llist));
              err|=sysMarkAlloc(sys,tlist->ptr,1);
              err|=sysMarkAlloc(sys,tlist->next,1);
            }
            break;

            case sysAllocTypePathlist:
            {
              PATHLIST* tpathlist;
              tpathlist=(PATHLIST*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(PATHLIST));
              err|=sysMarkAlloc(sys,tpathlist->path,1);
              err|=sysMarkAlloc(sys,tpathlist->next,1);
            }
            break;

            case sysAllocTypeNPath:
            {
              NPATH* tpath;
              tpath=(NPATH*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(NPATH));
              err|=sysMarkAlloc(sys,tpath->node,1);
              err|=sysMarkAlloc(sys,tpath->next,1);
            }
            break;

            case sysAllocTypeGPList:
            {
              GPLIST* tlist;
              tlist=(GPLIST*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(GPLIST));
              err|=sysMarkAlloc(sys,tlist->gp,1);
              err|=sysMarkAlloc(sys,tlist->next,1);
            }
            break;

            case sysAllocTypeSegList:
            {
              SEGLIST* tlist;
              unsigned int i, cnt;
              tlist=(SEGLIST*)ListPtr->AllocatedPtr;
              cnt=ListPtr->AllocatedSize/sizeof(SEGLIST);
              for (i=0; i<cnt; i++)
              {
                err|=sysMarkAlloc(sys,tlist->seg.segp,1);
                err|=sysMarkAlloc(sys,tlist->original,1);
                err|=sysMarkAlloc(sys,tlist->next,1);
                tlist++;
              }
            }
            break;

            case sysAllocTypePseudoSeg:
            {
              PSEUDO_SEG* tseg;
              tseg=(PSEUDO_SEG*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(PSEUDO_SEG));
              err|=sysMarkAlloc(sys,tseg->node[0],1);
              err|=sysMarkAlloc(sys,tseg->node[1],1);
              err|=sysMarkAlloc(sys,tseg->loops,1);
            }
            break;

            case sysAllocTypeExternal:
            {
              EXTERNAL* text;
              text=(EXTERNAL*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(EXTERNAL));
              err|=sysMarkAlloc(sys,text->source,1);
              err|=sysMarkAlloc(sys,text->indices,1);
              err|=sysMarkAlloc(sys,text->loops,1);
              err|=sysMarkAlloc(sys,text->name1,1);
              err|=sysMarkAlloc(sys,text->name2,1);
              err|=sysMarkAlloc(sys,text->portname,1);
              err|=sysMarkAlloc(sys,text->next,1);
            }
            break;

            case sysAllocTypeTree:
            {
              TREE* ttree;
              ttree=(TREE*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(TREE));
              err|=sysMarkAlloc(sys,ttree->loops,1);
              err|=sysMarkAlloc(sys,ttree->next,1);
            }
            break;

            case sysAllocTypeHoleList:
            {
              HoleList* thlist;
              thlist=(HoleList*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(HoleList));
              err|=sysMarkAlloc(sys,thlist->func,1);
              err|=sysMarkAlloc(sys,thlist->vals,1);
              err|=sysMarkAlloc(sys,thlist->next,1);
            }
            break;

            case sysAllocTypeOption:
            {
              Option* topt;
              topt=(Option*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(Option));
              err|=sysMarkAlloc(sys,topt->arg,1);
              err|=sysMarkAlloc(sys,topt->next,1);
            }
            break;

            case sysAllocTypeFilament:
            {
              FILAMENT* tfil;
              tfil=(FILAMENT*)ListPtr->AllocatedPtr;
              unsigned int cnt, i;
              cnt=ListPtr->AllocatedSize/sizeof(FILAMENT);
              for (i=0; i<cnt; i++)
              {
                err|=sysMarkAlloc(sys,tfil->segm,1);
                err|=sysMarkAlloc(sys,tfil->pchg,1);
                tfil++;
              }
            }
            break;

            case sysAllocTypeCharge:
            {
              charge* tchg;
              unsigned int cnt, i;
              cnt=ListPtr->AllocatedSize/sizeof(charge);
              tchg=(charge*)ListPtr->AllocatedPtr;
              for (i=0; i<cnt; i++)
              {
                err|=sysMarkAlloc(sys,tchg->next,1);
                err|=sysMarkAlloc(sys,tchg->surf,1);
                err|=sysMarkAlloc(sys,tchg->multipole,1);
                err|=sysMarkAlloc(sys,tchg->pos_dummy,1);
                err|=sysMarkAlloc(sys,tchg->neg_dummy,1);
                err|=sysMarkAlloc(sys,tchg->fil,1);
                tchg++;
              }
            }
            break;

            case sysAllocTypeSurface:
            {
              surface* tsurf;
              unsigned int cnt, i;
              tsurf=(surface*)ListPtr->AllocatedPtr;
              cnt=ListPtr->AllocatedSize/sizeof(surface);
              for(i=0; i<cnt; i++)
              {
                err|=sysMarkAlloc(sys,tsurf->title,1);
                err|=sysMarkAlloc(sys,tsurf->name,1);
                err|=sysMarkAlloc(sys,tsurf->group_name,1);
                err|=sysMarkAlloc(sys,tsurf->panels,1);
                err|=sysMarkAlloc(sys,tsurf->next,1);
                err|=sysMarkAlloc(sys,tsurf->prev,1);
                tsurf++;
              }
            }
            break;

            case sysAllocTypeNodeList:
            {
              NODELIST* tnlist;
              tnlist=(NODELIST*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(NODELIST));
              err|=sysMarkAlloc(sys,tnlist->name,1);
              err|=sysMarkAlloc(sys,tnlist->next,1);
            }
            break;


#ifdef SRW0814
            case sysAllocTypestEnt:
            {
              struct stEnt* h;
              h=(struct stEnt*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(struct stEnt))
              err|=sysMarkAlloc(sys,h->stNext,1);
              err|=sysMarkAlloc(sys,h->stTag,1);
              err|=sysMarkAlloc(sys,h->stData,1);
            }
            break;

            case sysAllocTypestTab:
            {
              struct stTab* h;
              h=(struct stTab*)ListPtr->AllocatedPtr;
              ASSERT (ListPtr->AllocatedSize==sizeof(struct stTab))
              err|=sysMarkAlloc(sys,h->entries,1);
            }
#endif

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
   printf("\n");
  {
    unsigned int tcnt;
    allocs=0;
    allocsize=0;
    ListPtr=sys->TopOfAllocationList;
    while (ListPtr!=NULL)
    {
      if ((ListPtr->Mark==0) && (ListPtr->AllocatedType!=0))
      {
        allocs++;
        allocsize+=ListPtr->AllocatedSize;
        printf("Orphant memory block typ %u @ %08x allocated in %s:%u ",ListPtr->AllocatedType, ListPtr->AllocatedPtr, ListPtr->AllocatedDebugFilename, ListPtr->AllocatedDebugLineNumber);
/*      if (ListPtr->AllocatedType==sysAllocTypeSpath)
        {
//            SPATH* tt = (SPATH*)ListPtr->AllocatedPtr;
//          printf("Source: %u\n", tt->src);
        }
*/
        printf("\n");
      }
      ListPtr=ListPtr->NextRecord;

    }
    if (allocs>0)
    {
      fprintf(stderr,"Orphant memory blocks %u, Size: %lu kB\n",allocs, allocsize>>10);
      exit (-1);
    }
  }
}

unsigned int sysMPIBlockNmbr;

/* */
/* Function to redirect data stream for load/save from/to FILE or MPI */
/* */
void sysComm(unsigned int mode, FILE* fop, void* ptr, unsigned long size)
{
  switch (mode)
  {
    #ifdef MPI
    case sysCommBcastMPI: // MPI Bcast
     MPI_Bcast(ptr, size/sizeof(char), MPI_CHAR, 0, *((MPI_Comm*)fop) );
     break;

    case sysCommtoMPI: // MPI Send
//     printf("MPI_Send sysMPIBlockNmbr %u size %u to %u \n", sysMPIBlockNmbr, size, *((int*)fop)); fflush (stdout);
     MPI_Send(ptr, size/sizeof(char), MPI_CHAR,
               *((int*)fop), sysMPIBlockNmbr++, MPI_COMM_WORLD);
     break;

    case sysCommfromMPI: // MPI Receive
//    printf("MPI_Recv sysMPIBlockNmbr %u size %u from %u \n", sysMPIBlockNmbr, size, *((int*)fop)); fflush (stdout);
     {
     MPI_Status status;
     MPI_Recv(ptr, size/sizeof(char), MPI_CHAR,
               *((int*)fop), sysMPIBlockNmbr++, MPI_COMM_WORLD, &status);
     }
     break;

    #endif

   case sysCommtoFile: // File save
     fwrite(ptr, size, 1, fop);
     break;

   case sysCommfromFile: // File read
     fread(ptr, size, 1, fop);
     break;

   default:
     fprintf(stderr, "spMatrixComm: Unknown communication target");
     exit (-1);
     break;
  }
 /**/
#ifdef MPI
  if (sysMPIBlockNmbr>=4095)
  {
    sysMPIBlockNmbr=0;
  }
#endif
}

/* Load, save, broadcast, send or receive matrix */
void sysCopyMatrix (double*** ptr, int rows, int cols, int size, int CommMode, FILE* fop)
{
    int MPIrank;
    int i;

    /* No pointer provided , so nothing to do */
    if (ptr==NULL)
    {
      return;
    }

   switch (CommMode)
   {
   case sysCommBcastMPI: // MPI Bcast
     #ifdef MPI
       ASSERT (fop!=NULL)
       MPI_Comm_rank ( *((MPI_Comm*)fop), &MPIrank);
     #else
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     break;

   case sysCommtoMPI:
     #ifndef MPI
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     ASSERT(fop!=NULL)
     ASSERT( (*((int*)fop))>=0)
     sysMPIBlockNmbr=1;
   case sysCommtoFile: // File save
     MPIrank=0;
     break;

   case sysCommfromMPI:
     #ifndef MPI
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     ASSERT(fop!=NULL)
     ASSERT( (*((int*)fop))>=0)
     sysMPIBlockNmbr=1;

   case sysCommfromFile: // File read
     MPIrank=1;
     break;

   default:
     /* Unknown communication target */
     ASSERT (1==0)
     break;
   }

   /* Copy link to first element */
   sysComm(CommMode, fop, (void*)ptr, sizeof(double**));


   /* Link list isn't there, nothing to do */
   if (*ptr==NULL)
   {
     return;
   }

   /* Allocate on destination side */
   if (MPIrank>0)
   {
     *ptr= (double **)MattAlloc(rows,sizeof(double *));
   }

   /* copy link list */
   sysComm(CommMode, fop, (void*)*ptr, rows*sizeof(double*));

   /* Copy content */
   for(i = 0; i < rows; i++)
   {

     if ((*ptr)[i]!=NULL)
     {
       /* Allocate on destination side */
       if (MPIrank>0)
       {
         (*ptr)[i] = (double *)MattAlloc(cols,size);
       }

       /* copy content */
       sysComm(CommMode, fop, (void*)(*ptr)[i], cols*size);
     }
   }
}


typedef struct memtable_t
{
	void* old;
	void* new;
	unsigned long length;
    unsigned int type;
} memtable_t;

#define sysUpdatePointer(ptr,memtable,memtableentries, hashtranstable) sysUpdatePointerFunc((void**)ptr,memtable, memtableentries, hashtranstable, __FILE__,__LINE__)
unsigned int sysUpdatePointerFunc(void** ptr, memtable_t* memtable, unsigned int memtableentries, unsigned long* hashtranstable, const char* filename, unsigned int linenumber)
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

  thash=sysCalcHash(*ptr);
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

  for (i=0;i<memtableentries; i++)
  {
    if ((lptr>=(unsigned long)memtable[i].old) && (lptr<(unsigned long)memtable[i].old+memtable[i].length))
    {
       lptr=lptr-(unsigned long)memtable[i].old+(unsigned long)memtable[i].new;
       *ptr=(void*)lptr;
       return 0;
    }
  }
  fprintf(stderr,"Pointer translation not successful %s:%u!\n",filename,linenumber);
  return 1;
*/
}

/* Load, save, broadcast, send or receive main data structure  */
void sysCopyindsys (SYS** psys, int CommMode, FILE* fop)
{
   int MPIrank;
   unsigned long i;

   ASSERT (psys!=NULL)

   switch (CommMode)
   {
   case sysCommBcastMPI: // MPI Bcast
     #ifdef MPI
       ASSERT (fop!=NULL)
       MPI_Comm_rank (*(MPI_Comm*)fop, &MPIrank);	/* get current process id */
     #else
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     break;

   case sysCommtoMPI:
     #ifndef MPI
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     sysMPIBlockNmbr=1;
   case sysCommtoFile: // File save
     MPIrank=0;
     break;

   case sysCommfromMPI:
     #ifndef MPI
       ASSERT(1==0)  /* Try to use MPI in a non MPI environment */
     #endif
     sysMPIBlockNmbr=1;

   case sysCommfromFile: // File read
     MPIrank=1;
     break;

   default:
     /* Unknown communication target */
     ASSERT (1==0)
     break;
   }

   if (MPIrank>0)
   {
     /* If there is an existing indsys structure -> release it */
     sysDestroy( psys );
   }

   /* Copy content */
   sysComm(CommMode, fop, (void*)psys, sizeof(SYS*));

   if ((*psys)==NULL)
   {
     /* Nothing to do */
     return;
   }

   {
      SYS* sys;
      sysAllocationListPtr  ListPtr;
      sysAllocationListHashPtr HashListPtr;
      memtable_t *memtable;
      unsigned long memtableentries, memtable_size;
      unsigned long hashtableentries[1<<sysCalcHashWidth],hashtableentries_all;
      unsigned long *hashtranstable,hashtranstable_size;

      if (MPIrank==0)
      {
        unsigned int cnt;

        sys = (*psys);
        ListPtr = sys->TopOfAllocationList;
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
        for (cnt=0; cnt< (1<<sysCalcHashWidth); cnt++ )
        {
           unsigned long hentries;
           sysAllocationListHashPtr HashListPtr;

           hentries=0;
           HashListPtr=sys->TopOfAllocationHashList[cnt];
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
        *psys = (SYS *)calloc(1, sizeof(SYS));
      }

      sysComm(CommMode, fop, (void*)&memtableentries, sizeof(unsigned long));
      sysComm(CommMode, fop, (void*)&hashtableentries_all, sizeof(unsigned long));

//      printf ("%u: MemTableEntries %lu HashTableEntries %lu\n", MPIrank, memtableentries, hashtableentries_all);
      memtable_size=(unsigned long)sizeof(struct memtable_t)*((unsigned long)memtableentries);
      memtable=malloc(memtable_size);
      ASSERT (memtable!=NULL);

      if (MPIrank==0)
      {
	    memtable_t *tmemtable;

        ListPtr = sys->TopOfAllocationList;
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

      hashtranstable_size=(unsigned long)sizeof(unsigned long) * (hashtableentries_all+(1u<<sysCalcHashWidth)+1);
      hashtranstable=malloc(hashtranstable_size);
      ASSERT (hashtranstable!=NULL)
      if (MPIrank==0)
      {
         unsigned int cnt;
         unsigned long entry;
         unsigned long *tlongptr;
         unsigned long tlcnt;
         entry=(1u<<sysCalcHashWidth)+1;
         tlongptr=hashtranstable;
         tlcnt=0;
         for (cnt=0; cnt< (1u<<sysCalcHashWidth); cnt++)
         {
            *(tlongptr++)=entry;
            entry+=hashtableentries[cnt];
            tlcnt++;
         }
         *(tlongptr++)=entry;
         tlcnt++;

         for (cnt=0; cnt< (1u<<sysCalcHashWidth); cnt++)
         {
           sysAllocationListHashPtr HashListPtr;

           // Assure consistency of table
           ASSERT(tlcnt==hashtranstable[cnt]);
           HashListPtr=sys->TopOfAllocationHashList[cnt];
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

      sysComm(CommMode, fop, (void*)memtable, memtable_size);
//      printf ("%u: MemTable copied\n", MPIrank); fflush(stdout);
      sysComm(CommMode, fop, (void*)hashtranstable, hashtranstable_size);
//      printf ("%u: HashTranstable copied (%u kB)\n", MPIrank,hashtranstable_size>>10); fflush(stdout);

      if (MPIrank>0)
      {
	    unsigned int memsize;
	    memtable_t *tmemtable;
        memsize=0;
        tmemtable=memtable;
        for (i=0; i<memtableentries; i++)
        {
           if ((tmemtable->type!=sysAllocTypeAllocationHashRecord) &&
               (tmemtable->type!=sysAllocTypeAllocationRecord))
           {
               ASSERT (tmemtable->length>0)

               /* Allocate memory */
               tmemtable->new=malloc(tmemtable->length);
               /* Register record and built immediatly local hash structure */
               sysRecordAllocationFunc(*psys,
                                       (char*)tmemtable->new,
                                        tmemtable->length,
                                        tmemtable->type,
                                        __FILE__,__LINE__);
               memsize+=tmemtable->length;
           }
           tmemtable++;
        }
//        printf ("%u: Memory allocated (%u kB)\n", MPIrank, memsize>>10);
      }
      {
    	memtable_t *tmemtable;
        tmemtable=memtable;
        for (i=0; i<memtableentries; i++)
        {
           #ifdef MPI
           if (CommMode==sysCommBcastMPI)
           {
             // Dirty: openmpi-1.7.3-1.fc20
             // For what ever reason more than 8192 Bcasts aren't allowed to be pending ...
             if ((i&2047)==0)
             {
               MPI_Barrier (*(MPI_Comm*)fop );
             }
           }
          #endif
          if ((tmemtable->type!=sysAllocTypeAllocationHashRecord) &&
               (tmemtable->type!=sysAllocTypeAllocationRecord))
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
             sysComm(CommMode, fop, (void*)temp, tmemtable->length);

             if (MPIrank>0)
             {
             err=0;
	         /* Update pointers in copied content */
             switch (tmemtable->type)
             {
             case sysAllocTypeGeneric:
               /* Nothing to do */
               break;

             case sysAllocTypePtrArray:
             {
               unsigned int tcount;
               int** tptr;
               unsigned int i;
               tcount=tmemtable->length/sizeof(int*);
               tptr=(int**)temp;
               for (i=0; i<tcount; i++)
               {
                 err|=sysUpdatePointer(&(*tptr),memtable,memtableentries,hashtranstable);
                 tptr++;
               }
             }
             break;

             case AllocTypessystem:
             {
               ssystem *tsys;
               tsys=(ssystem*)temp;
               ASSERT(tmemtable->length==sizeof(ssystem))
               err|=sysUpdatePointer(&tsys->cond_names,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->q,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->p,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->panels,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->cubes,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->multilist,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->locallist,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->directlist,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->precondlist,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->revprecondlist,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->is_dummy,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->is_dielec,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsys->indsys,memtable,memtableentries,hashtranstable);
             }
             break;

             case AllocTypecube:
             {
               cube *tcube;
               tcube=(cube*)temp;
               ASSERT(tmemtable->length==sizeof(cube))
               err|=sysUpdatePointer(&tcube->mnext,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->upnumeles,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->upvects,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->multi,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->upmats,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->is_dummy,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->is_dielec,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->lnext,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->downnumeles,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->downvects,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->local,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->downmats,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->interList,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->enext,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->evalnumeles,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->evalvects,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->eval,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->evalmats,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->eval_isQ2P,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->dnext,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->pnext,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->rpnext,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->directnumeles,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->directq,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->directmats,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->precondmats,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->directlu,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->precond,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->prevectq,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->prevectp,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->nbr_is_dummy,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->chgs,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->nbrs,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->kids,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcube->parent,memtable,memtableentries,hashtranstable);
             }
             break;

             case AllocTypecharge:
             {
               charge *tcharge;
               tcharge=(charge*)temp;
               ASSERT(tmemtable->length==sizeof(charge))
               err|=sysUpdatePointer(&tcharge->next,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcharge->multipole,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcharge->surf,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcharge->pos_dummy,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcharge->neg_dummy,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tcharge->fil,memtable,memtableentries,hashtranstable);
             }
             break;

             case AllocTypesurface:
             {
               surface *tsurface;
               tsurface=(surface*)temp;
               ASSERT(tmemtable->length==sizeof(surface))
               err|=sysUpdatePointer(&tsurface->title,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsurface->name,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsurface->panels,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsurface->group_name,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsurface->next,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tsurface->prev,memtable,memtableentries,hashtranstable);
             }
             break;

             case AllocTypeName:
             {
               Name *tname;
               tname=(Name*)temp;
               ASSERT(tmemtable->length==sizeof(Name))
               err|=sysUpdatePointer(&tname->name,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tname->next,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tname->alias_list,memtable,memtableentries,hashtranstable);
             }
             break;


            case sysAllocTypeindsys:
            { SYS *tsys;
              unsigned int i;
              tsys=(SYS*)temp;
              ASSERT (tmemtable->length==sizeof(SYS))

              err|=sysUpdatePointer(&tsys->planes,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->endplane,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->segment,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->endseg,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->nodes,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->endnode,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->pseudo_nodes,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->meshsect,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->diagL,memtable,memtableentries,hashtranstable);
              tsys->sparMatrix=NULL;
              err|=sysUpdatePointer(&tsys->Alist,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->M,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->Mlist,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->m_info,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->Mtrans,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->Precond,memtable,memtableentries,hashtranstable);;
              err|=sysUpdatePointer(&tsys->trees,memtable,memtableentries,hashtranstable);;
              err|=sysUpdatePointer(&tsys->externals,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->Z,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->R,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->MtZM,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->FinalY,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->resids,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->resid_real,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->resid_imag,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->niters,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->opts,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->Ar,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->Br,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->Cr,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->Dr,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->title,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->chglist,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->sys,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->EmptyMElement,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->EmptySEGLIST,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->EmptyFILAMENT,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->EmptyCHARGE,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->EmptySURFACE,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->EmptySEGMENT,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->EmptyNODE,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tsys->EmptySPATH,memtable,memtableentries,hashtranstable);
#ifdef SRW0814
              err|=sysUpdatePointer(&tsys->tab.entries,memtable,memtableentries,hashtranstable);
#endif
              /* Copy memory management links */
              for (i=0; i<(1u<<sysCalcHashWidth); i++)
              {
                tsys->TopOfAllocationHashList[i]=(*psys)->TopOfAllocationHashList[i];
              }
              tsys->TopOfAllocationList=(*psys)->TopOfAllocationList;
              tsys->EmptyAllocHashRecord=(*psys)->EmptyAllocHashRecord;
              tsys->EmptyAllocRecord=(*psys)->EmptyAllocRecord;

              free (*psys);
              *psys=tsys;
            }
            break;

            case sysAllocTypeindopts:
            {
             ind_opts *tind_opts;
             tind_opts=(ind_opts*)temp;
             ASSERT (tmemtable->length==sizeof(ind_opts))
             err|=sysUpdatePointer(&tind_opts->portlist,memtable,memtableentries,hashtranstable);
             err|=sysUpdatePointer(&tind_opts->suffix,memtable,memtableentries,hashtranstable);
             err|=sysUpdatePointer(&tind_opts->fname,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypestrlist:
            {
             strlist *tstrlist;
             tstrlist=(strlist*)temp;
             ASSERT (tmemtable->length==sizeof(strlist))
             err|=sysUpdatePointer(&tstrlist->str,memtable,memtableentries,hashtranstable);
             err|=sysUpdatePointer(&tstrlist->next,memtable,memtableentries,hashtranstable);
            }
            case sysAllocTypeNode:
            {
              NODES* tnode;
              unsigned int i, cnt;
              tnode=(NODES*)temp;
              cnt=tmemtable->length/sizeof(NODES);
              for (i=0;i<cnt; i++)
              {
                err|=sysUpdatePointer(&tnode->name,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tnode->equiv,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tnode->to_end,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tnode->connected_segs,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tnode->gp,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tnode->treeptr,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tnode->gp_node,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tnode->next,memtable,memtableentries,hashtranstable);
                tnode++;
              }
            }
            break;

            case sysAllocTypeSegment:
            {
              SEGMENT* tseg;
              unsigned int i, cnt;
              tseg=(SEGMENT*)temp;
              cnt=tmemtable->length/sizeof(SEGMENT);
              for (i=0; i<cnt; i++)
              {
                err|=sysUpdatePointer(&tseg->name,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tseg->widthdir,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tseg->node[0],memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tseg->node[1],memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tseg->filaments,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tseg->loops,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tseg->next,memtable,memtableentries,hashtranstable);
                tseg++;
              }
            }
            break;

            case sysAllocTypeGroundplane:
            {
              GROUNDPLANE *tgp;
              tgp=(GROUNDPLANE*)temp;
              ASSERT (tmemtable->length==sizeof(GROUNDPLANE))
              err|=sysUpdatePointer(&tgp->name,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->grid1,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->grid2,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->segs1,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->segs2,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->pnodes,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->lastnode,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->innode,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->outnode,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->usernodes,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->fake_seg_list,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->list_of_holes,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->list_of_contacts,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->usernode_coords,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->filename,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->nonuni,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->indsys,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypePseudoNode:
            {
              PSEUDO_NODE* tnode;
              tnode=(PSEUDO_NODE*)temp;
              ASSERT (tmemtable->length==sizeof(PSEUDO_NODE))
              err|=sysUpdatePointer(&tnode->node,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tnode->name,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tnode->next,memtable,memtableentries,hashtranstable);
            }
            break;


            case sysAllocTypeSpath:
            {
              SPATH* tspath;
              unsigned int i, cnt;
              tspath=(SPATH*)temp;
              cnt=tmemtable->length/sizeof(SPATH);
              for (i=0; i<cnt; i++)
              {
                err|=sysUpdatePointer(&tspath->seg.segp,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tspath->next,memtable,memtableentries,hashtranstable);
                tspath++;
              }
            }
            break;

            case sysAllocTypeMElement:
            {
              MELEMENT* tmelement;
              unsigned int tcnt, tcount;
              tmelement=(MELEMENT*)temp;
              tcount=tmemtable->length/sizeof(MELEMENT);
              for (tcnt=0; tcnt<tcount; tcnt++)
              {
                err|=sysUpdatePointer(&tmelement->fil,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tmelement->mnext,memtable,memtableentries,hashtranstable);
                tmelement++;
              }
            }
            break;

            case sysAllocTypeContactList:
            {
              ContactList* tcontact;
              tcontact=(ContactList*)temp;
              ASSERT (tmemtable->length==sizeof(ContactList))
              err|=sysUpdatePointer(&tcontact->func,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tcontact->name,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tcontact->vals,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tcontact->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeBi:
            {
              Bi* tbi;
              tbi=(Bi*)temp;
              ASSERT (tmemtable->length==sizeof(Bi))
              err|=sysUpdatePointer(&tbi->child1,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tbi->child2,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeGrid_2d:
            {
              Grid_2d* tgrid_2d;
              tgrid_2d=(Grid_2d*)temp;
              ASSERT (tmemtable->length==sizeof(Grid_2d))
              err|=sysUpdatePointer(&tgrid_2d->kids,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeInd_list:
            {
              int_list* tlist;
              tlist=(int_list*)temp;
              ASSERT (tmemtable->length==sizeof(int_list))
              err|=sysUpdatePointer(&tlist->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeNonuni_gp:
            {
              Nonuni_gp *tgp;
              tgp=(Nonuni_gp*)temp;
              ASSERT (tmemtable->length==sizeof(Nonuni_gp));
              err|=sysUpdatePointer(&tgp->root_cell,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->z_pts,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->thick,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->z_c,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->nodelist,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgp->grndp,memtable,memtableentries,hashtranstable);
            }
            break;


            case sysAllocTypeG_edges:
            {
              G_edges *tedge;
              unsigned int i;
              tedge=(G_edges*)temp;
              ASSERT (tmemtable->length==sizeof(G_edges))
              for (i=0; i<NUM_E_CELLS; i++)
              {
                err|=sysUpdatePointer(&tedge->cells[i],memtable,memtableentries,hashtranstable);
              }
              err|=sysUpdatePointer(&tedge->children,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeGcell:
            {
              Gcell *tgcell;
              unsigned int i,cnt,j ;
              tgcell=(Gcell*)temp;
              cnt=tmemtable->length/sizeof(Gcell);
              for (j=0; j<cnt; j++)
              {
                err|=sysUpdatePointer(&tgcell->children,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tgcell->parent,memtable,memtableentries,hashtranstable);
                for (i=0; i<MAX(NUMEDGES,NUMNODES); i++)
                {
                  err|=sysUpdatePointer(&tgcell->bndry.edges[i],memtable,memtableentries,hashtranstable);
                }
                tgcell++;
              }
            }
            break;

            case sysAllocTypeG_nodes:
            {
              G_nodes *tgnode;
              unsigned int i;
              tgnode=(G_nodes*)temp;
              ASSERT (tmemtable->length==sizeof(G_nodes));
              for (i=0; i<NUM_N_CELLS; i++)
              {
                err|=sysUpdatePointer(&tgnode->cells[i],memtable,memtableentries,hashtranstable);
              }
              for (i=0; i<NUMADJ; i++)
              {
                err|=sysUpdatePointer(&tgnode->adjacent[i],memtable,memtableentries,hashtranstable);
              }
              err|=sysUpdatePointer(&tgnode->e_segs,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgnode->n_segs,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgnode->prev,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tgnode->next,memtable,memtableentries,hashtranstable);
            }
            break;


            case sysAllocTypeLlist:
            {
              Llist *tlist;
              tlist=(Llist*)temp;
              ASSERT (tmemtable->length==sizeof(Llist))
              err|=sysUpdatePointer(&tlist->ptr,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tlist->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypePathlist:
            {
              PATHLIST* tpathlist;
              tpathlist=(PATHLIST*)temp;
              ASSERT (tmemtable->length==sizeof(PATHLIST))
              err|=sysUpdatePointer(&tpathlist->path,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tpathlist->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeNPath:
            {
              NPATH* tpath;
              tpath=(NPATH*)temp;
              ASSERT (tmemtable->length==sizeof(NPATH))
              err|=sysUpdatePointer(&tpath->node,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tpath->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeGPList:
            {
              GPLIST* tlist;
              tlist=(GPLIST*)temp;
              ASSERT (tmemtable->length==sizeof(GPLIST))
              err|=sysUpdatePointer(&tlist->gp,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tlist->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeSegList:
            {
              SEGLIST* tlist;
              unsigned int i, cnt;
              tlist=(SEGLIST*)temp;
              cnt=tmemtable->length/sizeof(SEGLIST);
              for (i=0; i<cnt; i++)
              {
               err|=sysUpdatePointer(&tlist->seg.segp,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tlist->original,memtable,memtableentries,hashtranstable);
               err|=sysUpdatePointer(&tlist->next,memtable,memtableentries,hashtranstable);
               tlist++;
              }
            }
            break;

            case sysAllocTypePseudoSeg:
            {
              PSEUDO_SEG* tseg;
              tseg=(PSEUDO_SEG*)temp;
              ASSERT (tmemtable->length==sizeof(PSEUDO_SEG))
              err|=sysUpdatePointer(&tseg->node[0],memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tseg->node[1],memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tseg->loops,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeExternal:
            {
              EXTERNAL* text;
              text=(EXTERNAL*)temp;
              ASSERT (tmemtable->length==sizeof(EXTERNAL))
              err|=sysUpdatePointer(&text->source,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&text->indices,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&text->loops,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&text->name1,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&text->name2,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&text->portname,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&text->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeTree:
            {
              TREE* ttree;
              ttree=(TREE*)temp;
              ASSERT (tmemtable->length==sizeof(TREE))
              err|=sysUpdatePointer(&ttree->loops,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&ttree->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeHoleList:
            {
              HoleList* thlist;
              thlist=(HoleList*)temp;
              ASSERT (tmemtable->length==sizeof(HoleList))
              err|=sysUpdatePointer(&thlist->func,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&thlist->vals,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&thlist->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeOption:
            {
              Option* topt;
              topt=(Option*)temp;
              ASSERT (tmemtable->length==sizeof(Option))
              err|=sysUpdatePointer(&topt->arg,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&topt->next,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypeFilament:
            {
              FILAMENT* tfil;
              unsigned int i,cnt;
              cnt=tmemtable->length/sizeof(FILAMENT);
              tfil=(FILAMENT*)temp;
              for (i=0; i<cnt; i++)
              {
                err|=sysUpdatePointer(&tfil->segm,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tfil->pchg,memtable,memtableentries,hashtranstable);
                tfil++;
              }
            }
            break;

            case sysAllocTypeCharge:
            {
              charge* tchg;
              unsigned int cnt, i;
              cnt=tmemtable->length/sizeof(charge);
              tchg=(charge*)temp;
              for (i=0; i<cnt; i++)
              {
                err|=sysUpdatePointer(&tchg->next,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tchg->surf,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tchg->multipole,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tchg->pos_dummy,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tchg->neg_dummy,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tchg->fil,memtable,memtableentries,hashtranstable);
                tchg++;
              }
            }
            break;

            case sysAllocTypeSurface:
            {
              surface* tsurf;
              unsigned int cnt, i;
              tsurf=(surface*)temp;
              cnt=tmemtable->length/sizeof(surface);
              for (i=0; i<cnt; i++)
              {
                err|=sysUpdatePointer(&tsurf->title,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tsurf->name,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tsurf->group_name,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tsurf->panels,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tsurf->next,memtable,memtableentries,hashtranstable);
                err|=sysUpdatePointer(&tsurf->prev,memtable,memtableentries,hashtranstable);
                tsurf++;
              }
            }
            break;

            case sysAllocTypeNodeList:
            {
              NODELIST* tnlist;
              tnlist=(NODELIST*)temp;
              ASSERT (tmemtable->length==sizeof(NODELIST))
              err|=sysUpdatePointer(&tnlist->name,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&tnlist->next,memtable,memtableentries,hashtranstable);
            }
            break;

#ifdef SRW0814
            case sysAllocTypestEnt:
            {
              struct stEnt* h;
              h=(struct stEnt*)temp;
              ASSERT (tmemtable->length==sizeof(struct stEnt))
              err|=sysUpdatePointer(&h->stNext,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&h->stTag,memtable,memtableentries,hashtranstable);
              err|=sysUpdatePointer(&h->stData,memtable,memtableentries,hashtranstable);
            }
            break;

            case sysAllocTypestTab:
            {
              struct stTab* h;
              h=(struct stTab*)temp;
              ASSERT (tmemtable->length==sizeof(struct stTab))
              err|=sysUpdatePointer(&h->entries,memtable,memtableentries,hashtranstable);
            }
#endif // SRW0814

             default:
               fprintf(stderr,"CopyindsysoverMPI: Unknown AllocatedType %u\n",tmemtable->type);
               ASSERT (1==0);
             }
             if  (err>0)
             {
               fprintf(stderr,"Block typ %u @ %08x allocated contains links outside of memory allocation table -> Links corrupt!\n", tmemtable->type, tmemtable->old);
               ASSERT (1==0);
             }
            }
          }
          tmemtable++;
        }
      }
      //printf ("%u: MemTable linked content copied \n", MPIrank);

      free(memtable);
      free(hashtranstable);
   }

}


unsigned int sysCalcHash(void* ptr)
{

  unsigned long addr;
  unsigned int thash;

  addr=((unsigned long)ptr)>>sysCalcHashGridSize;

  thash = ((addr                      )&((1u<<sysCalcHashWidth)-1));
  thash^= ((addr>>(  sysCalcHashWidth))&((1u<<sysCalcHashWidth)-1));
  thash^= ((addr>>(2*sysCalcHashWidth))&((1u<<sysCalcHashWidth)-1));
  thash^= ((addr>>(3*sysCalcHashWidth))&((1u<<sysCalcHashWidth)-1));
  thash^= ((addr>>(4*sysCalcHashWidth))&((1u<<sysCalcHashWidth)-1));
  thash^= ((addr>>(5*sysCalcHashWidth))&((1u<<sysCalcHashWidth)-1));
//  thash^= ((addr>>(6*sysCalcHashWidth))&((1u<<sysCalcHashWidth)-1));

  return thash;
}

void sysAnalyze2 (SYS* indsys)
{
  sysAllocationListHashPtr tptr;


  unsigned int hashtabpop[(1u<<sysCalcHashWidth)];
  unsigned int i,records,cnt;
  unsigned int population[100];


  if (indsys==NULL)
   return;

  for (i=0; i<100; i++)
  {
    population[i]=0;
  }


  records=0; cnt=0;
  for (i=0; i<(1u<<sysCalcHashWidth); i++)
  {
    hashtabpop[i]=0;
    tptr=indsys->TopOfAllocationHashList[i];
    while (tptr!=NULL)
    {
      ASSERT (tptr->hashtable==i)
      ASSERT (tptr->AllocRecord!=NULL)
      ASSERT (tptr->AllocRecord->AllocatedPtr!=NULL)
      ASSERT (tptr->AllocRecord->AllocatedSize!=0)
      cnt++;
      hashtabpop[i]++;
      if (tptr->PrevPeerRecord==NULL)
      {
       records++;
       population[tptr->AllocRecord->AllocatedType]++;
      }
      tptr=tptr->NextRecord;
    }
  }
/*
  for (i=0; i<100; i++)
  if (population[i]>0)
  printf("Type %u: %u\n",i,population[i]);
*/
/*
    {
    unsigned int cnt1,i,j;
    cnt1=0;
    for (i=0; i<MAX( ((1u<<sysCalcHashWidth)>>4),1); i++)
    {
     for (j=0; j<MIN(16, 1u<<sysCalcHashWidth); j++)
     {
       printf("%03u ", hashtabpop[cnt1++]);
     }
     printf("\n");
    }
    printf("%u Records %u Hashs\n",records, cnt);
   }
*/
//  if ((cnt>>(2*9) )>0)
  {
    unsigned int i;
    i=4;
    while ((cnt>>(2*i))>0)
     i++;
    if (i>sysCalcHashWidth)
    {
     printf("Hashwidth %u recommended\nUpdate sysCalcHashWidth in indsys.h and recompile\n",i);
    }
  }

/*
  printf ("EmptySPATH %08x\n",indsys->EmptySPATH);
  {
    SPATH* temp;
    char* a, *b;
    temp=indsys->EmptySPATH;

    a=(char*)&temp->next;
    b=(char*)temp;

    printf("Start %08x %08x %u %u\n",temp, &temp->next,a-b,sizeof(char));
  }
*/

}
