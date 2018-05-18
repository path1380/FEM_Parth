  program main
#include <petsc/finclude/petsc.h>
  use type_defs
  use petsc
  use petscvec
  implicit none
  
  PetscInt n
  PetscErrorCode ierr
  PetscBool flg
  PetscScalar one,two,three,dot
  PetscReal norm,rdot
  Vec x,y,w
  PetscOptions options

  n = 20
  one = 1.0
  two = 2.0
  three = 3.0

  call MPI_Init(ierr)
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  !call PetscOptionsCreate(options,ierr)
  !call VecCreate(PETSC_COMM_SELF,x,ierr)

  write(*,*) "Program runs!"
!contains

end program main
