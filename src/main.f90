  program main
#include <petsc/finclude/petsc.h>
  use type_defs
  use petsc
  use petscvec
  implicit none
  
  PetscInt n,sz
  PetscErrorCode ierr
  PetscBool flg
  PetscScalar one,two,three,dot,a22
  PetscReal norm,rdot
  Vec x,y,w,b
  Mat A
  PetscOptions options

  n = 2
  one = 1.0
  a22 = -1.0
  two = 2.0
  three = 3.0

  call MPI_Init(ierr)
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call PetscOptionsCreate(options,ierr)
  call PetscOptionsGetInt(options,PETSC_NULL_CHARACTER,'-n',n,flg,ierr)
  call PetscOptionsDestroy(options,ierr)

  call VecCreate(PETSC_COMM_SELF,b,ierr)
  call VecSetSizes(b,PETSC_DECIDE,n,ierr)
  call VecSetFromOptions(b,ierr)

  call MatCreate(PETSC_COMM_SELF,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
  call MatSetFromOptions(A,ierr)
  call MatSetup(A,ierr)

  !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
  call MatSetValues(A,1,1,1,1,a22,INSERT_VALUES,ierr)
  call MatSetValues(A,1,0,1,0,one,INSERT_VALUES,ierr)
  call MatSetValues(A,1,1,1,0,one,INSERT_VALUES,ierr)
  call MatSetValues(A,1,0,1,1,one,INSERT_VALUES,ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  call MatView(A,PETSC_VIEWER_STDOUT_SELF,ierr)

  call MatScale(A,two,ierr)

  call MatView(A,PETSC_VIEWER_STDOUT_SELF,ierr)

  !call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call VecGetSize(x,sz,ierr)
  !call PetscIntView(1,sz,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call VecDuplicate(b,y,ierr)
  !call VecDuplicate(b,w,ierr)

  !call VecSet(x,one,ierr)
  !call VecSet(y,two,ierr)

  !call VecDot(x,y,dot,ierr)
  !rdot = PetscRealPart(dot)

  !write(6,*) rdot

  !call VecScale(x,two,ierr)

  !call VecNorm(x,NORM_2,norm,ierr)

  !call VecAXPY(y,three,x,ierr)

  !call VecView(y,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !write(*,*) norm
  !call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call PetscScalarView(1,dot,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)

  write(*,*) "Program runs!"
!contains

end program main
