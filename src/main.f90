  program main
#include <petsc/finclude/petsc.h>
  use type_defs
  use petsc
  use petscvec
  implicit none
  
  PetscInt n,sz,i,j
  PetscErrorCode ierr
  PetscBool flg
  PetscScalar zero,one,two,three,dot,a22
  PetscReal norm,rdot
  Vec x,y,w,b
  Mat A
  KSP ksp
  PetscOptions options

  n = 3
  zero = 0.0
  one = 1.0
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

  call VecSetValues(b,1,0,three,INSERT_VALUES,ierr)
  call VecSetValues(b,1,1,two,INSERT_VALUES,ierr)
  call VecSetValues(b,1,2,one,INSERT_VALUES,ierr)

  !call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)

  call MatCreate(PETSC_COMM_SELF,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
  call MatSetFromOptions(A,ierr)
  call MatSetup(A,ierr)

  !call MatView(A,PETSC_VIEWER_STDOUT_SELF,ierr)

  do i=0,n-1
     do j=0,n-1
        if (i <= j) then
           call MatSetValues(A,1,i,1,j,one,INSERT_VALUES,ierr)
        else
           call MatSetValues(A,1,i,1,j,zero,INSERT_VALUES,ierr)
        end if
     end do
  end do

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  !call MatView(A,PETSC_VIEWER_STDOUT_SELF,ierr)
  
  call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
  call KSPSetFromOptions(ksp,ierr)
  
  call KSPSetOperators(ksp,A,A,ierr)

  call KSPSetType(ksp,KSPCG,ierr)

  call KSPSolve(ksp,b,x,ierr)

  call VecView(x,PETSC_VIEWER_STDOUT_SELF,ierr)

  !call KSPDestroy(ksp,ierr)
  !call VecDestroy(b,ierr)
  !call VecDestroy(x,ierr)
  !call MatDestroy(A,ierr)
  write(*,*) "Program runs!"
!contains

end program main
