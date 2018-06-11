module solver_context
#include <petsc/finclude/petsc.h>
  use type_defs
  use petsc
  implicit none
  TYPE MatCtx
     real(kind=dp), dimension(4,4) :: local_matrix
  END TYPE MatCtx

end module solver_context
