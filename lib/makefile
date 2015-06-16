# set the C compiler, for C++ use g++

# use this setting on most sane machines
CC = gcc

# use this setting on the engin elf machines
#CC = gcc296

# set the path to the include directory
INCDIR =../include

# set the flags for the C and C++ compiler to give lots of warnings
CFLAGS = -I$(INCDIR) -O2 -Wall -Wstrict-prototypes -Wnested-externs -Wmissing-prototypes -Wmissing-declarations
CPPFLAGS = $(CFLAGS)

# library tool defs
AR = ar
ARFLAGS = rv
RANLIB = ranlib

# path to the object file directory
ODIR = obj

# path to the lib directory
LIBDIR =.

# path to the bin directory
BINDIR =../bin

# put all of the relevant include files here
_DEPS = graphics.h

# convert them to point to the right place
DEPS = $(patsubst %,$(INCDIR)/%,$(_DEPS))

# put a list of all the object files (with .o endings)
_COMMON = color.o image.o ppmIO.o blackhole.o rk4.o

# convert them to point to the right place
COMMON = $(patsubst %,$(ODIR)/%,$(_COMMON))


# patterns for compiling source code
$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.C $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# output library
OUTPUT_LIB = libimageIO.a

# make rules
all: $(OUTPUT_LIB)

$(OUTPUT_LIB) : $(COMMON) $(DEPS)
	rm -f $@
	$(AR) $(ARFLAGS) $@ $(COMMON)
	$(RANLIB) $@

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 

