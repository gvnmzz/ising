#
SHELL=/bin/bash
SOURCE=./
PREFIX=./

OBJS := main.o
MODULES_OBJS := mtwist.o

# Provide any options and flags for the compiler here
FFLAGS= -O2 -Wall -I$(PREFIX)

# This is the compiler
F = gfortran
FC = $(F) $(FFLAGS)
LD = $(F)

# Set name of final executable file here
EXENAME = ising2d

default:$(EXENAME)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# compile and generate dependency info
%.o: %.f*
	$(FC) -c $< -o $@

# link to generate executable file
$(EXENAME): $(MODULES_OBJS) $(OBJS)
	$(LD) $(MODULES_OBJS) $(OBJS) -o $(EXENAME)


# remove compilation products
.PHONY : clean

clean:
	rm -f $(EXENAME) *.o *.d *.mod *~ *.pdf play.pbs.* *.res
