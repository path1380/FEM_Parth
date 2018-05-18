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
  KSP ksp
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

  call VecDuplicate(b,x,ierr)

  call VecSet(b,one,ierr)

  !call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)

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

  !call MatView(A,PETSC_VIEWER_STDOUT_SELF,ierr)
  
  call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
  call KSPSetFromOptions(ksp,ierr)
  
  call KSPSetOperators(ksp,A,A,ierr)

  call KSPSolve(ksp,b,x,ierr)

  call VecView(x,PETSC_VIEWER_STDOUT_SELF,ierr)

  call KSPDestroy(ksp,ierr)
  call VecDestroy(b,ierr)
  call VecDestroy(x,ierr)
  call MatDestroy(A,ierr)
  write(*,*) "Program runs!"
!contains

end program main
