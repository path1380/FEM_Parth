module solver_context
#include <petsc/finclude/petsc.h>
  use type_defs
  use basic_structures
  use petsc
  implicit none
  TYPE MatCtx
     real(kind=dp), dimension(4,4) :: local_matrix
     type(problem_data) :: prob
     type(numerics_data) :: num
  END TYPE MatCtx

end module solver_context
