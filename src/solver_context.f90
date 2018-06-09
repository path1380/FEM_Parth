module solver_context
#include <petsc/finclude/petsc.h>
  use type_defs
  use petsc
  implicit none
  TYPE MatCtx
     PetscReal :: lambda,kappa
     PetscReal :: h
  END TYPE MatCtx

end module solver_context
