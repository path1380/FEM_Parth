  program main
#include <petsc/finclude/petsc.h>
  use type_defs
  use petsc
  use petscvec
  implicit none
  
  PetscInt n,sz
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
  call PetscOptionsCreate(options,ierr)
  call PetscOptionsGetInt(options,PETSC_NULL_CHARACTER,'-n',n,flg,ierr)
  call PetscOptionsDestroy(options,ierr)

  call VecCreate(PETSC_COMM_SELF,x,ierr)
  call VecSetSizes(x,PETSC_DECIDE,n,ierr)
  call VecSetFromOptions(x,ierr)
  call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call VecGetSize(x,sz,ierr)

  write(*,*) "Program runs!"
!contains

end program main
