# Compiler and flags
Exe = ../Kernel_Check
FC = nvfortran
FCFLAGS = -O3
LDFLAGS = -fopenmp -Mfree -Mpreprocess 
#LIBS = -lcufft



# Test mode flag (set TEST=1 to enable)
ifdef TEST
    FCFLAGS += -DTEST
    Exe := $(Exe)-TEST
endif


Files = global.f90 Kernel.f90 Kernel_Check.f90
Obs = $(Files:.f90=.o)

# Linking

                                                                     
$(Exe): $(Obs)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(Exe) $(Obs) 

#$(LIBS)

# Compilation
%.o: %.f90
	$(FC) $(FCFLAGS) $(LDFLAGS) -c $<


clean:
	rm -f $(Obs) $(Exe) *.mod

               
.PHONY: clean