# Makefile for FEM_Parth

PETSC_DIR ?= /home/parth/petsc
PETSC_ARCH ?= arch-linux2-c-debug

SRC = src
BIN = bin
FOPT = -O0 -g
F90FLAGS = -fbounds-check -Wall -fbacktrace $(FOPT)
FFLAGS = -cpp -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
_SRCS = type_defs.f90 basic_structures.f90 Input_function.f90 basis_function.f90 numbering_convention_defn.f90 legendre_module.f90 solver_context.f90 solver_context_interfaces.f90 local_mat.f90 global_mat.f90 rs_to_xy.f90 approx_fun.f90 main.f90
SRCS = $(patsubst %,$(SRC)/%,$(_SRCS))
_OBJS = $(_SRCS:.f90=.o)
OBJECTS = $(patsubst %,$(BIN)/%,$(_OBJS))
EXECUTABLE = test.x
# EXECUTABLE = test_pis.x

.PHONY: clean compile build run

include $(PETSC_DIR)/lib/petsc/conf/variables

all : build

$(OBJECTS) : $(BIN)/%.o : $(SRC)/%.f90
	$(FC) $(FFLAGS) $(F90FLAGS) -c -o $@ $(FC_MODULE_OUTPUT_FLAG)$(@D) $^

$(EXECUTABLE) : $(OBJECTS)
	$(FC_LINKER) -o $@ $^ -llapack $(PETSC_LIB)

compile: $(OBJECTS)

build: $(EXECUTABLE)

run: $(EXECUTABLE)
	$(MPIEXEC) -np 1 ./$(EXECUTABLE) >> output.txt

clean::
	rm -f $(OBJECTS) $(EXECUTABLE) $(BIN)/*.mod
	rm -f *.o *.mod
