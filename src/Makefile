# Start of customization part.

# Installation directories.
PREFIX = /usr/local
INCDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib

# Libraries to build.
LIBS = libclbfgsb3.a libclbfgsb3.so

# C and FORTRAN compilers.
CC = gcc
FC = gfortran

# Compiler flags.
CFLAGS = -Wall -O3 -mavx2 -mfma -ffast-math -fPIC
#CFLAGS = -Wall -O2 -fPIC
FFLAGS = $(CFLAGS)

# Linker flags.
LDFLAGS = -lm

# Flags to build a shared library.
SHLIB_FLAGS = -shared

# Commands.
RM = rm -f
COPY = cp -a
MAKEDIR = mkdir -p

# End of customization part.
#------------------------------------------------------------------------------

# Source directory.
srcdir = .

# L-BFGS-B code.
LBFGSB_SRCDIR = $(srcdir)/../lbfgsb-3.0
LBFGSB_SRCS = \
    $(LBFGSB_SRCDIR)/blas.f \
    $(LBFGSB_SRCDIR)/lbfgsb.f \
    $(LBFGSB_SRCDIR)/linpack.f \
    $(LBFGSB_SRCDIR)/timer.f

OBJS = \
    blas.o \
    clbfgsb.o \
    lbfgsb.o \
    linpack.o \
    timer.o

TESTS = \
    clbfgsb_test1 \
    clbfgsb_test2 \
    clbfgsb_test3

TEST_OUTPUTS = \
    clbfgsb_test1.out \
    clbfgsb_test2.out \
    clbfgsb_test3.out

default: $(LIBS) $(TESTS) $(TEST_OUTPUTS)

install: $(LIBS)
	$(MAKEDIR) "$(INCDIR)"
	$(COPY) lbfgsb.h "$(INCDIR)"
	$(MAKEDIR) "$(LIBDIR)"
	for src in $(LIBS); do $(COPY) $$src "$(LIBDIR)"; done

clean:
	$(RM) *.o *~

dist-clean: clean
	$(RM) $(LIBS) $(TESTS) $(TEST_OUTPUTS) iterate.dat

check: $(TEST_OUTPUTS)

libclbfgsb3.a: $(OBJS)
	ar rv $@ $^

libclbfgsb3.so: $(OBJS)
	$(FC) $(SHLIB_FLAGS) -o $@ $^ $(LDFLAGS)

%.out: %
	./$< | sed -e 's/\([0-9]\)[eE]\([-+][0-9]\)/\1D\2/g' >$@

clbfgsb_test1: clbfgsb_test1.o $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

clbfgsb_test1.o: $(srcdir)/clbfgsb_test1.c $(srcdir)/lbfgsb.h
	$(CC) -I$(srcdir) $(CFLAGS) -o $@ -c $<

clbfgsb_test2: clbfgsb_test2.o $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

clbfgsb_test2.o: $(srcdir)/clbfgsb_test2.c $(srcdir)/lbfgsb.h
	$(CC) -I$(srcdir) $(CFLAGS) -o $@ -c $<

clbfgsb_test3: clbfgsb_test3.o $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)
clbfgsb_test3.o: $(srcdir)/clbfgsb_test3.c $(srcdir)/lbfgsb.h
	$(CC) -I$(srcdir) $(CFLAGS) -o $@ -c $<

clbfgsb.o: $(srcdir)/clbfgsb.c $(srcdir)/lbfgsb.h
	$(CC) -I$(srcdir) $(CFLAGS) -o $@ -c $<

blas.o: $(LBFGSB_SRCDIR)/blas.f
	$(FC) $(FFLAGS) -o $@ -c $<

lbfgsb.o: $(LBFGSB_SRCDIR)/lbfgsb.f
	$(FC) $(FFLAGS) -o $@ -c $<

linpack.o: $(LBFGSB_SRCDIR)/linpack.f
	$(FC) $(FFLAGS) -o $@ -c $<

timer.o: $(LBFGSB_SRCDIR)/timer.f
	$(FC) $(FFLAGS) -o $@ -c $<

.PHONY: clean dist-clean check default install
