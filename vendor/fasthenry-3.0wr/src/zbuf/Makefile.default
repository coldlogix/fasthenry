# Makefile for
#  zbuf - capacitance or inductance visualization picture

# These CFLAGS should be fine for anything other than MinGW.
CFLAGS = -O -g -DFOUR

# For MinGW, there is no getrusage, and we'll skip a gethostname call
# to avoid having to link with winsock.
#CFLAGS = -O -g -DFOUR -DNO_RUSAGE -DNO_GETHOSTNAME

# We're going to assume gcc, you'll need to comment this if a non-compatible
# compiler is used.
CC = gcc

# assumed places for unaliased executables - must be same in all Makefile.*
SHELL = /bin/sh
AWK = awk
RM = /bin/rm

MAKETYPE = other
# normally overidden by master Makefile
TYPE = other

INCL = .
IO = .
DIR = .
MUL = .
ALC = .
COM = .
GEN = .
BIN = ../../bin
TOP = ..

OBJS = $(MUL)/mulGlobal.o $(DIR)/extras.o \
$(MUL)/zbuf.o \
$(ALC)/uglieralloc.o \
$(IO)/patran.o $(IO)/quickif.o \
$(IO)/input.o \
$(IO)/zbufInOut.o $(IO)/zbuf2fastcap.o $(IO)/zbufProj.o $(IO)/zbufSort.o \
$(IO)/dump_struct.o $(IO)/savemat_mod.o

SRCS = $(MUL)/mulGlobal.c $(DIR)/extras.c \
$(MUL)/zbuf.c \
$(ALC)/uglieralloc.c \
$(IO)/patran.c $(IO)/quickif.c \
$(IO)/input.c \
$(IO)/zbufInOut.c $(IO)/zbuf2fastcap.c $(IO)/zbufProj.c $(IO)/zbufSort.c \
$(IO)/dump_struct.c $(IO)/savemat_mod.c

LIBS = -lm

zbuf: $(OBJS)
	$(CC) -o zbuf $(CFLAGS) $(OBJS) $(LIBS)
	mv zbuf $(BIN)/zbuf

clean:
	$(RM) -f *.o

$(OBJS) : $(INCL)/mulStruct.h $(INCL)/mulGlobal.h $(INCL)/patran.h \
$(INCL)/resusage.h Makefile
$(SRCS) : $(INCL)/mulStruct.h $(INCL)/mulGlobal.h $(INCL)/patran.h \
$(INCL)/resusage.h Makefile
$(IO)/quickif.c $(IO)/quickif.o : $(IO)/quickif.h Makefile
