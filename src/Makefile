# Compiler and flags
Exe = ../TBSBIEM-v1.1.0
FC = nvfortran
FCFLAGS = -acc -cuda -gpu=managed -fast -O3
 
LDFLAGS = -fopenmp -Mfree -Mpreprocess
LIBS = -lcufft


# Test mode flag (set TEST=1 to enable)
ifdef TEST
    FCFLAGS += -DTEST
    Exe := $(Exe)-TEST
endif

ifdef Minfo
    FCFLAGS  += -g -Minfo=accel
endif

Files = global.cuf Kernel.cuf cuFFT2D.cuf norms.cuf main.cuf 
Obs = $(Files:.cuf=.o)

# Linking

                                                                     
$(Exe): $(Obs)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(Exe) $(Obs) $(LIBS)

# Compilation
%.o: %.cuf
	$(FC) $(FCFLAGS) $(LDFLAGS) -c $<

              
clean:
	rm -f $(Obs) $(Exe) *.mod

               
.PHONY: clean