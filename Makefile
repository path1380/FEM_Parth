# Makefile for FEM_Parth

PETSC_DIR ?= /home/parth/petsc
PETSC_ARCH ?= arch-linux2-c-debug

MPIEXEC = mpiexec
FC = mpif90
LD = $(FC)
SRC = src
BIN = bin
F90FLAGS = -fbounds-check -Wall -fbacktrace -g
FFLAGS = -O3 -cpp -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
_SRCS = type_defs.f90 basic_structures.f90 Input_function.f90 basis_function.f90 numbering_convention_defn.f90 legendre_module.f90 local_mat.f90 global_mat.f90 rs_to_xy.f90 approx_fun.f90 solver_context.f90 solver_context_interfaces.f90 main.f90
SRCS = $(patsubst %,$(SRC)/%,$(_SRCS))
_OBJS = $(_SRCS:.f90=.o)
OBJECTS = $(patsubst %,$(BIN)/%,$(_OBJS))
EXECUTABLE = test.x
# EXECUTABLE = test_pis.x

.PHONY: clean compile build run

include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules


compile:
	$(FC) $(FFLAGS) -c $(SRCS) -g
	mv *.mod $(BIN)
	mv *.o $(BIN)

build: $(OBJECTS)
	$(LD) -o $(EXECUTABLE) $(OBJECTS) -llapack $(PETSC_LIB)

run: $(EXECUTABLE)
	$(MPIEXEC) -np 1 ./$(EXECUTABLE) >> output.txt

clean::
	rm -f $(OBJECTS) $(EXECUTABLE) $(BIN)/*.mod
	rm -f *.o *.mod
