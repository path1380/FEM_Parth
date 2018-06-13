module solver_context_interfaces
#include <petsc/finclude/petsc.h>
  use type_defs
  use petsc
  use solver_context

  ! ----------------------------------------------------
         INTERFACE MatCreateShell
           SUBROUTINE MatCreateShell(comm,mloc,nloc,m,n,ctx,mat,ierr)
             USE solver_context
             MPI_Comm :: comm
             PetscInt :: mloc,nloc,m,n
             TYPE(MatCtx) :: ctx
             Mat :: mat
             PetscErrorCode :: ierr
           END SUBROUTINE MatCreateShell
         END INTERFACE MatCreateShell
  ! ----------------------------------------------------

  ! ----------------------------------------------------
         INTERFACE MatShellSetContext
           SUBROUTINE MatShellSetContext(mat,ctx,ierr)
             USE solver_context
             Mat :: mat
             TYPE(MatCtx) :: ctx
             PetscErrorCode :: ierr
           END SUBROUTINE MatShellSetContext
         END INTERFACE MatShellSetContext
  ! ----------------------------------------------------

  ! ----------------------------------------------------
         INTERFACE MatShellGetContext
           SUBROUTINE MatShellGetContext(mat,ctx,ierr)
             USE solver_context
             Mat :: mat
             TYPE(MatCtx),  POINTER :: ctx
             PetscErrorCode :: ierr
           END SUBROUTINE MatShellGetContext
         END INTERFACE MatShellGetContext

! ----------------------------------------------------
         INTERFACE PCShellSetContext
           SUBROUTINE PCShellSetContext(pc,ctx,ierr)
             USE solver_context
             PC :: pc
             TYPE(MatCtx) :: ctx
             PetscErrorCode :: ierr
           END SUBROUTINE PCShellSetContext
         END INTERFACE PCShellSetContext
  ! ----------------------------------------------------

  ! ----------------------------------------------------
         INTERFACE PCShellGetContext
           SUBROUTINE PCShellGetContext(pc,ctx,ierr)
             USE solver_context
             PC :: pc
             TYPE(MatCtx),  POINTER :: ctx
             PetscErrorCode :: ierr
           END SUBROUTINE PCShellGetContext
         END INTERFACE PCShellGetContext

end module solver_context_interfaces
