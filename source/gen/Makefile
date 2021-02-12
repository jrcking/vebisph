#==================================================
# 
# Dr. Xiaohu Guo
# xiaohu.guo@stfc.ac.uk
#
#==================================================
#need to specify architecture unfortunately
arch_bit := sixtyfour


SHELL=/bin/sh
MAKE=gmake
FC=gfortran

ifeq ($(arch_bit), thirtytwo)
FCFLAGS=-W -O3 -m32 -g -fbounds-check
FCA = -m32
else ifeq ($(arch_bit), sixtyfour)
FCFLAGS=-W -O3 -m64 -g -fbounds-check
FCA = -m64
endif

VPATH =../para  

OBJS= kind_parameters.o common_parameter.o global_variables.o datclass.o

default: gen2D

gen2D: $(OBJS)
	$(FC) $(FCA) -o $@ $(OBJS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod gen2D isphgen
