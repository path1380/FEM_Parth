!================================================================================================================
!
! File: main.f90
! Brief description: Basic solver for weak form of \int_{\Omega} v (u - f) d \vec{x} = 0 for all v. Bilinear
!                    basis function used. The code can be extended to higher order basis functions easily.
!
! Detailed description: 
!
!
! Author: Parth Thakkar
!
!================================================================================================================

program main
#include <petsc/finclude/petsc.h>
  use type_defs
  use basic_structures
  use basis_function
  use Input_function
  use numbering_convention_defn
  use legendre_module
  use local_mat
  use global_mat
  use rs_to_xy
  use approx_fun
  use petsc
  implicit none

  real(kind=dp), allocatable, dimension(:) :: qnodes,qweights
  real(kind=dp), dimension(4,4) :: A_local
  real(kind=dp), dimension(4) :: b_local,b_local1,b_local2,b_local3,b_local4
  real(kind=dp) :: test_x,test_y
  
  real(kind=dp), allocatable, dimension(:,:) :: A_global  
  real(kind=dp), allocatable, dimension(:) :: b_global,b_test,f_val

  integer :: num_divs_x,num_divs_y,num_nodes
  
  type(node) :: test_node
  type(element) :: test_element,test_element1,test_element2,test_element3,test_element4

  integer :: i,nx,ny
  real(kind=dp) :: f_approx,xi,yi,hx,hy,length,width

!MPI variables

  integer :: nprocs,myid

!LAPACK variables

  integer, allocatable, dimension(:) :: IPIV
  integer :: INFO

!Petsc variables

  integer, allocatable, dimension(:) :: row_ind,col_ind
  PetscInt n
  PetscErrorCode ierr
  PetscBool flg
  PetscOffset i_b,i_soln

  Vec b,soln
  Mat A

  KSP ksp
  PC pc

  PetscOptions options

  PetscViewer viewer
  PetscDraw draw

!=================Petsc Initializing=============================================
  !Initializing Petsc
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  !Giving number of processors
  call MPI_COMM_SIZE(PETSC_COMM_WORLD,nprocs,ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD,myid,ierr)

  call PetscOptionsCreate(options,ierr)
  call PetscOptionsGetInt(options,PETSC_NULL_CHARACTER,'-n',num_nodes,flg,ierr)
  !call PetscOptionsDestroy(options,ierr)
!=================Petsc Initializing=============================================
 

  !Giving problem data input
  call Input_problem_data(prob_data_test)

  length = prob_data_test%domain_length
  width = prob_data_test%domain_width

  !Giving numerics data input
  call Input_numerics_data(num_data_test)

  num_divs_x = num_data_test%num_divs_x
  num_divs_y = num_data_test%num_divs_y
  num_nodes = (num_divs_x+1)*(num_divs_y+1)

  n = num_nodes

!=====================Petsc declarations========================================
  
  allocate(row_ind(0:num_nodes-1))
  allocate(col_ind(0:num_nodes-1))

  do i=0,num_nodes-1
     row_ind(i) = i
     col_ind(i) = i
  end do

  !Creating Petsc Vectors b and soln
  call VecCreate(PETSC_COMM_WORLD,b,ierr)
  call VecSetSizes(b,PETSC_DECIDE,n,ierr)
  call VecSetFromOptions(b,ierr)

  call VecDuplicate(b,soln,ierr)

  !Creating Petsc Matrix
  call MatCreate(PETSC_COMM_WORLD,A,ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
  call MatSetFromOptions(A,ierr)
  call MatSetup(A,ierr)

  !Creating Petsc Viewer
  !call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER,'Hello world!',&
  !     &PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,viewer,ierr)

!=====================Petsc declarations========================================

  hx = length/num_divs_x
  hy = width/num_divs_y

  !Allocating Global matrices
  allocate(A_global(num_nodes,num_nodes))

  allocate(b_global(num_nodes))
  allocate(qnodes(0:num_data_test%num_quadrature_nodes-1))
  allocate(qweights(0:num_data_test%num_quadrature_nodes-1))

  allocate(b_test(num_nodes))

  allocate(f_val(num_nodes))
  

  allocate(IPIV(num_nodes))

  !test_node = solve_node(1, prob_data_test, num_data)
  !test_element1 = solve_element(6, prob_data_test, num_data_test)
  !write(*,*) r(0.25_dp,test_element1)
  !write(*,*) s(0.8_dp,test_element1)
  !test_element2 = solve_element(2, prob_data_test, num_data_test)
  !test_element3 = solve_element(3, prob_data_test, num_data_test)
  !test_element4 = solve_element(4, prob_data_test, num_data_test)

  !write(*,*) test_node%neigh_elt_nos
  !write(*,*) test_element%nodes%coord(1)

  !call build_local_A(prob_data_test,num_data_test, A_local)

  call build_global_matrices(A_global, b_global, prob_data_test, num_data_test)

  !write(*,*) nint
  !write(*,*) A_global - transpose(A_global)
  !write(*,*) b_global
  !write(*,*) A_local - transpose(A_local)
  !write(*,*) b_local
  !do i=1,121
  !   write(*,*) b_global(i)
  !end do  

  !===========================Solution using LAPACK============================

  !call DGETRF(num_nodes,num_nodes,A_global,num_nodes,IPIV,INFO)
  !call DGETRS('N',num_nodes,1,A_global,num_nodes,IPIV,b_global,num_nodes,INFO)

  !f_approx = approx_eval(0.2_dp,0.2_dp,b_global,prob_data_test,num_data_test)

  !===========================Solution using LAPACK============================


  !write(*,*) f_approx
  
! CODE USED TO CHECK ERROR AT POINTS OTHER THAN NODES
  !do i=1,num_nodes
  !   b_test(i) = 0.0_dp
  !   ny = num_divs_y+1-mod(i-1,num_divs_y+1)
  !   nx = 1 + ((i-1)/(num_divs_y+1))
  !   xi = hx*(nx-1)
  !   yi = hy*(ny-1)
  !   b_test(i) = approx_eval(xi,yi,b_global,prob_data_test,num_data_test)
  !end do
! CODE USED TO CHECK ERROR AT POINTS OTHER THAN NODES

!==============Using Petsc-ksp for solving=============
  !write(*,*) A_global

  !call VecCreateSeqWithArray(PETSC_COMM_SELF,PETSC_DECIDE,num_nodes,b_global,b,ierr)
  !call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
  !write(*,*) b_global

  call VecSetValues(b,num_nodes,col_ind,b_global,INSERT_VALUES,ierr)
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

  call MatSetValues(A,num_nodes,row_ind,num_nodes,col_ind,transpose(A_global),INSERT_VALUES,ierr)
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  !call MatView(A,PETSC_VIEWER_STDOUT_SELF,ierr)
  
  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  call KSPSetFromOptions(ksp,ierr)
  
  call KSPSetOperators(ksp,A,A,ierr)

  call KSPSetType(ksp,KSPCG,ierr)

  call KSPGetPC(ksp,pc,ierr)
  call PCSetType(pc,PCJACOBI,ierr)

  !do
    call KSPSolve(ksp,b,soln,ierr)
  !end do

  !call KSPSetFromOptions(ksp,ierr)
  

  !call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
  !write(*,*) b_global
!==============Using Petsc-ksp for solving=============

!======================Plotting using Petsc Viewer===============================

  !call PetscViewerDrawGetDraw(viewer,0,draw,ierr)

!======================Plotting using Petsc Viewer===============================

  do i=1,num_nodes
     f_val(i) = 0.0_dp
     ny = num_divs_y+1-mod(i-1,num_divs_y+1)
     nx = 1 + ((i-1)/(num_divs_y+1))
     xi = hx*(nx-1)
     yi = hy*(ny-1)
     f_val(i) = f(xi,yi)
  end do
  !write(*,*) log(dble(num_divs_x)),log(maxval(abs(b_global - f_val)))

!==============Destroying and Finalizing Petsc objects==

  call KSPDestroy(ksp,ierr)
  call MatDestroy(A,ierr)
  call VecDestroy(b,ierr)
  call VecDestroy(soln,ierr)
  call PetscOptionsDestroy(options,ierr)

  call PetscFinalize(ierr)

!==============Destroying and Finalizing Petsc objects==
end program main
