# Makefile for sonic method
#
# set up variables that we will use.   The first COMPILER variable decides which compiler/linker 
# and associated flags will be chosen.
COMPILER = gnu
#COMPILER = intel

ifeq ($(COMPILER),intel) 
  CXX = icc 
  CFLAGS = -O3 -xSSSE3 
  LINK = icc 
  LDFLAGS = 
endif

ifeq ($(COMPILER),gnu) 
  CXX = /usr/bin/gcc
  CFLAGS = -O3 -mssse3 -DANSI
  LINK = $(CXX) 
  LDFLAGS = -lm
endif

HEADERS = ss.h
SRCS = main.c ss.c
OBJECTS = $(SRCS:.C=.o)
DISTFILES = $(SRCS) $(HEADERS) Makefile


# here are the targets.  The target is listed before the colon, and
# after the colon are the files that it depends on.  make will ensure
# that the dependencies are updated before executing the command(s)
# on the following lines.

# this target tells make how to LINK the test_2d_array.o object file with
# the necessary libraries and startup code (even though most complications are masked by
# our use of the compiler to do the job of linking).
sonic: $(OBJECTS)
	$(LINK) $(LDFLAGS) -o $@ $^ 

# we explicitly tell Make that the "clean" and "dist" targets never really have a file 
# associated with them. 
.PHONY: clean dist

# here is a (phony) target to clean up our source directory
clean:
	-rm -f *.o
	-rm -f sonic

# here is a (phony) target to create a distribution tar file from our source code. 
dist:
	tar cvfz dist.tar.gz $(DISTFILES)


# these implicit suffix rules tell Make how to create .o files from .C or .c files.  Basically,
# this gives the command that we will use to COMPILE.
.C.o:
	$(CXX) $(CFLAGS) -c $<
.c.o:
	$(CXX) $(CFLAGS) -c $<
