# Makefile
FC      := gfortran
FFLAGS  := -ffixed-line-length-none -std=legacy -fPIC -Wl,--no-relax
SRC     := $(wildcard *.f)
OBJ     := $(SRC:.f=.o)

.PHONY: all debug release clean

all: release

release: FFLAGS += -O3
release: qseis_stress.bin

debug:   FFLAGS += -g -O0
debug:   qseis_stress_dbg.bin

qseis_stress.bin: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

qseis_stress_dbg.bin: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o qseis_stress.bin qseis_stress_dbg.bin
