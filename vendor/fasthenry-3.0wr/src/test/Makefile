CC=gcc
CCFLAGSADD=-DMPI -I/usr/include/openmpi-x86_64 -m64 -g
LFLAGSADD=-L/usr/lib64/openmpi/lib -lmpi -lmpi_cxx

#CCFLAGSADD=-DMPI -I/usr/include/mpich-x86_64 -Wimplicit-function-declaration
#LFLAGSADD=-lmpich -L/usr/lib64/mpich/lib
CFLAGS = -O -DFOUR  $(CCFLAGSADD)
# -Wuninitialized -Wimplicit-int
SHELL = /bin/sh
RM = /bin/rm
BIN = ../../bin

SPLIB = ../fasthenry/sparse/sparse.a

OBJS = test.o

MOBJS =

NONUNIOBJS =

HEADER = cmplx.h
MULHEAD =
NONUNIHEAD =

test:	$(OBJS) $(MOBJS) $(NONUNIOBJS) $(SPLIB)
	$(CC) -o test $(CFLAGS) $(OBJS) $(MOBJS) $(NONUNIOBJS) $(SPLIB) -lm $(LFLAGSADD)
	mv test $(BIN)/test

$(SPLIB):
	cd ../fasthenry/sparse; make CC=$(CC) CFLAGSADD="$(CCFLAGSADD)" LFLAGSADD="$(LFLAGSADD)"

clean:
	$(RM) -f *.o

$(OBJS):	$(HEADER) $(MULHEAD)

$(MOBJS):	$(MULHEAD)

$(NONUNIOBJS):  $(HEADER) $(MULHEAD) $(NONUNIHEAD)
