
# determine availibility of MPI and set compiler and linker options accordingly
ifneq ($(shell which gcc 2>/dev/null),)
CC=gcc
endif

ifeq ($(shell which mpicc 2>/dev/null),)
CCFLAGDMPI =
else
CC=mpicc
CCFLAGDMPI = -DMPI
endif

CCFLAGSADD=$(CCFLAGDMPI)  -g -Wimplicit-function-declaration
LFLAGSADD =

CFLAGS = -O -DFOUR  $(CCFLAGSADD)
# -Wuninitialized -Wimplicit-int
SHELL = /bin/sh
RM = /bin/rm
BIN = ../../bin

SPLIB = sparse/sparse.a

OBJS = induct.o gmres.o savemat_mod.o readGeom.o joelself.o writefastcap.o \
       SetupMulti.o calcp.o SetupComputePsi.o mulSetup.o BreakupSeg.o \
       Precond.o addgroundplane.o findpaths.o fillM.o cx_ludecomp.o \
       parse_command_line.o default_opts.o Prec_cost.o hole.o \
       dist_betw_fils.o mutual.o newPrecond.o deg_mutual.o \
       barnoldi.o regurgitate.o memmgmt.o

MOBJS = mulGlobal.o mulDo.o mulMulti.o mulLocal.o mulMats.o direct.o \
        uglieralloc.o capsolve.o

NONUNIOBJS = find_nonuni_path.o read_tree.o contact.o

HEADER = induct.h cmplx.h resusage.h
MULHEAD = mulStruct.h mulGlobal.h patran.h resusage.h
NONUNIHEAD = gp.h

fasthenry:	$(OBJS) $(MOBJS) $(NONUNIOBJS) $(SPLIB)
	$(CC) -o fasthenry $(CFLAGS) $(OBJS) $(MOBJS) $(NONUNIOBJS) $(SPLIB) -lm $(LFLAGSADD)
	mv fasthenry $(BIN)/fasthenry

$(SPLIB):
	cd sparse; make CC=$(CC) CFLAGSADD="$(CCFLAGSADD)" LFLAGSADD="$(LFLAGSADD)"

clean:
	$(RM) -f *.o
	cd sparse; make clean

$(OBJS):	$(HEADER) $(MULHEAD)

$(MOBJS):	$(MULHEAD)

$(NONUNIOBJS):  $(HEADER) $(MULHEAD) $(NONUNIHEAD)
