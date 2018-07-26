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
  use solver_context_interfaces
  implicit none

  real(kind=dp), allocatable, dimension(:) :: qnodes,qweights
  real(kind=dp), dimension(4,4) :: A_local
  !real(kind=dp), dimension(4) :: b_local1,b_local2,b_local3,b_local4
  !real(kind=dp) :: test_x,test_y
  !real(kind=dp), dimension(4) :: b_local_test
  !real(kind=dp), dimension(4) :: local_soln
  integer, dimension(4) :: loc_col_ind
  !real(kind=dp), dimension(4) :: op,ret
  !integer, dimension(4) :: elt_glonos
  
  !real(kind=dp), allocatable, dimension(:,:) :: A_global  
  real(kind=dp), allocatable, dimension(:) :: b_global,f_val,b_global_shell
  !real(kind=dp), allocatable, dimension(:) :: b_test

  integer :: num_divs_x,num_divs_y,num_nodes,num_elements
  
  !type(node) :: test_node
  !type(element) :: test_element,test_element1,test_element2,test_element3,test_element4

  integer :: i,nx,ny
  !integer :: j
  integer :: big_loop_variable
  real(kind=dp) :: xi,yi,hx,hy,length,width
  !real(kind=dp) :: f_approx

!MPI variables

  integer :: nprocs,myid

!LAPACK variables

  !integer, allocatable, dimension(:) :: IPIV
  !integer :: INFO

!Petsc variables

  integer, allocatable, dimension(:) :: row_ind,col_ind
  integer :: n_iter_newton
  PetscInt n
  PetscErrorCode ierr
  !PetscBool flg
  !PetscOffset i_b,i_soln
  PetscScalar tol
  !PetscScalar norm_local_delta_u
  PetscScalar norm_delta_u
  !PetscScalar temp_norm

  TYPE(MatCtx) :: ctxA
  !TYPE(MatCtx),POINTER :: ctxA_pt,ctxA1

  Vec b,soln_iter,b_newton,soln_init,delta_u,soln_prev,temp_vec
  !Vec local_delta_u,local_soln_iter,local_soln_prev,local_b_newton
  !Vec temp_op_arg,temp_ret_val
  !Vec soln
  !Vec b_local
  Mat A
  Mat A_global_shell
  Mat A_local_shell
  !Mat A_local_petsc
  !Mat pshellmat

  !KSP ksp,ksp_iter
  KSP ksp_iter_shell
  !KSP ksp_local_shell
  !PC pc
  PC pc_shell
  !PC pc_shell_local

  !type(element) :: test_element
  !type(element) :: elt


  !PetscViewer viewer
  !PetscDraw draw

!=================Petsc Initializing=============================================
  !Initializing Petsc
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  !Giving number of processors
  call MPI_COMM_SIZE(PETSC_COMM_WORLD,nprocs,ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD,myid,ierr)

  !call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',num_nodes,flg,ierr)
!=================Petsc Initializing=============================================
 

  !Giving problem data input
  call Input_problem_data(prob_data_test)

  length = prob_data_test%domain_length
  width = prob_data_test%domain_width

  !Giving numerics data input
  call Input_numerics_data(num_data_test)

!===========Testing build_local_A_petsc==================================

  !call MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,4,4,&
  !                   &PETSC_NULL_INTEGER,A_local_petsc,ierr)

  !call MatCreate(PETSC_COMM_WORLD,A_local_petsc,ierr)
  !call MatSetSizes(A_local_petsc,PETSC_DECIDE,PETSC_DECIDE,4,4,ierr)
  !call MatSetFromOptions(A_local_petsc,ierr)
  !call MatSetType(A_local_petsc,MATSHELL,ierr)
  !call MatSetup(A_local_petsc,ierr)

  !call build_local_A_petsc(prob_data_test,num_data_test, A_local_petsc)

  !call MatAssemblyBegin(A_local_petsc,MAT_FINAL_ASSEMBLY,ierr)
  !call MatAssemblyEnd(A_local_petsc,MAT_FINAL_ASSEMBLY,ierr)

  !call MatView(A_local_petsc,PETSC_VIEWER_STDOUT_WORLD,ierr)

  
  !call MatShellGetContext(A_global_shell,ctxA_pt,ierr)

  !write(*,*) ctxA_pt%local_matrix

  !do i=1,4
  !   write(*,*) A_local(i,:)
  !end do

!  stop 123

!===========Testing build_local_A_petsc==================================

!====================MAIN CODE FOR SOLVING=====================================

!BIG LOOP FOR POST-PROCESSING

  big_loop_variable = 2

!  do while(big_loop_variable <= 550)

     num_data_test%num_divs_x = big_loop_variable
     num_data_test%num_divs_y = big_loop_variable

     num_divs_x = num_data_test%num_divs_x
     num_divs_y = num_data_test%num_divs_y

     num_nodes = (num_divs_x+1)*(num_divs_y+1)
     num_elements = num_divs_x*num_divs_y

     n = num_nodes

     call build_local_A(prob_data_test,num_data_test, A_local)

     !do i=1,4
     !   write(*,*) A_local(i,:)
     !end do


     ctxA%local_matrix = A_local
     ctxA%prob = prob_data_test
     ctxA%num = num_data_test


     call MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,&
       &num_nodes,num_nodes,ctxA,A_global_shell,ierr)

     call MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,&
       &4,4,ctxA,A_local_shell,ierr)


     call MatCreate(PETSC_COMM_WORLD,A_global_shell,ierr)
     call MatSetSizes(A_global_shell,PETSC_DECIDE,PETSC_DECIDE,num_nodes,num_nodes,ierr)
     call MatSetFromOptions(A_global_shell,ierr)
     call MatSetType(A_global_shell,MATSHELL,ierr)

     call MatSetup(A_global_shell,ierr)
     call MatSetup(A_local_shell,ierr)

     call MatShellSetContext(A_global_shell,ctxA,ierr)
     call MatShellSetContext(A_local_shell,ctxA,ierr)

     call MatShellSetOperation(A_global_shell,MATOP_MULT,MyMult,ierr)
     call MatShellSetOperation(A_local_shell,MATOP_MULT,MyMult_local,ierr)

!=====================Petsc declarations========================================
  
     allocate(row_ind(0:num_nodes-1))
     allocate(col_ind(0:num_nodes-1))

!=========Testing the matrix multiplication operation definition=============

!     num_divs_x = num_data_test%num_divs_x
!     num_divs_y = num_data_test%num_divs_y

!     num_nodes = (num_divs_x+1)*(num_divs_y+1)

     !allocate(row_ind(0:num_nodes-1))


!     call VecCreate(PETSC_COMM_WORLD,temp_op_arg,ierr)
!     call VecSetSizes(temp_op_arg,PETSC_DECIDE,num_nodes,ierr)
!     call VecSetFromOptions(temp_op_arg,ierr)

!     call VecDuplicate(temp_op_arg,temp_ret_val,ierr)

!     do i=0,num_nodes-1
  
!        call VecSetValue(temp_op_arg,i,1.0_dp,INSERT_VALUES,ierr)

!        call VecAssemblyBegin(temp_op_arg,ierr)
!        call VecAssemblyEnd(temp_op_arg,ierr)

!        call MatMult(A_global_shell,temp_op_arg,temp_ret_val,ierr)

!        call VecView(temp_ret_val,PETSC_VIEWER_STDOUT_WORLD,ierr)

!        call VecSet(temp_op_arg,0.0_dp,ierr)
!        call VecSet(temp_ret_val,0.0_dp,ierr)

!     end do

!     allocate(A_global(num_nodes,num_nodes))
!     allocate(b_global(num_nodes))

!     call build_global_matrices(A_global, b_global, prob_data_test, num_data_test)

!     do i=1,num_nodes
!        write(*,*) A_global(i,:)
!     end do

!=========Testing the matrix multiplication operation definition=============

!====================Testing MyMult_local================================

     !call VecCreate(PETSC_COMM_WORLD,temp_op_arg,ierr)
     !call VecSetSizes(temp_op_arg,PETSC_DECIDE,4,ierr)
     !call VecSetFromOptions(temp_op_arg,ierr)

     !call VecDuplicate(temp_op_arg,temp_ret_val,ierr)

     !do i=0,3
  
     !   call VecSetValue(temp_op_arg,i,1.0_dp,INSERT_VALUES,ierr)

     !   call VecAssemblyBegin(temp_op_arg,ierr)
     !   call VecAssemblyEnd(temp_op_arg,ierr)

     !   call MatMult(A_local_shell,temp_op_arg,temp_ret_val,ierr)

     !   call VecView(temp_ret_val,PETSC_VIEWER_STDOUT_WORLD,ierr)

     !   call VecSet(temp_op_arg,0.0_dp,ierr)
     !   call VecSet(temp_ret_val,0.0_dp,ierr)

     !end do

     !do i=1,4
     !   write(*,*) A_local(i,:)
     !end do


!====================Testing MyMult_local================================

!====================Testing build_local_b_vec===========================

     !do i=1,4
     !  test_element = solve_element(i,prob_data_test,num_data_test)
     !  call build_local_b(prob_data_test,num_data_test,test_element,b_local_test)
     !  write(*,*) b_local_test

     !  call VecCreate(PETSC_COMM_WORLD,b_local,ierr)
     !  call VecSetSizes(b_local,PETSC_DECIDE,4,ierr)
     !  call VecSetFromOptions(b_local,ierr)
     !  call build_local_b_vec(prob_data_test,num_data_test,test_element,b_local)

     !  call VecView(b_local,PETSC_VIEWER_STDOUT_SELF,ierr)
     !end do

!====================Testing build_local_b_vec===========================

!====================Testing local preconditioner========================

     !call VecCreate(PETSC_COMM_WORLD,temp_op_arg,ierr)
     !call VecSetSizes(temp_op_arg,PETSC_DECIDE,4,ierr)
     !call VecSetFromOptions(temp_op_arg,ierr)

     !call VecSet(temp_op_arg,1.0_dp,ierr)

     !call VecDuplicate(temp_op_arg,temp_ret_val,ierr)

     !call PCCreate(PETSC_COMM_WORLD,pc_shell_local,ierr)
     !call PCSetType(pc_shell_local,PCSHELL,ierr)
     !call PCShellSetContext(pc_shell_local,ctxA,ierr)
     !call PCShellSetApply(pc_shell,PC_Shell_Jacobi,ierr)

     !call PC_Shell_Jacobi_local(pc_shell_local,temp_op_arg,temp_ret_val,ierr)
     
     !call VecView(temp_ret_val,PETSC_VIEWER_STDOUT_SELF,ierr)

     

!====================Testing local preconditioner========================

!====================Testing libCEED_Mult================================

     !do i=1,4
     !   op(i_) = 0.0_dp
     !end do

     !do i=1,4
     !   op(i) = 1.0_dp
     !   call libCEED_Mult(prob_data_test,num_data_test,op,ret)
     !   write(*,*) ret
     !   write(*,*) A_local(:,i)
     !   do j=1,4
     !      op(j) = 0.0_dp
     !   end do
     !end do

!====================Testing libCEED_Mult================================


     do i=0,num_nodes-1
        row_ind(i) = i
        col_ind(i) = i
     end do

     do i=1,4
        loc_col_ind(i) = i-1
     end do

  !Creating Petsc Vectors b and soln
     call VecCreate(PETSC_COMM_WORLD,b,ierr)
     call VecSetSizes(b,PETSC_DECIDE,n,ierr)
     call VecSetFromOptions(b,ierr)

  !   call VecDuplicate(b,soln,ierr)

  !Creating vectors for iterative solution
     call VecDuplicate(b,soln_iter,ierr)
     call VecDuplicate(b,b_newton,ierr)
     call VecDuplicate(b,soln_init,ierr)
     call VecDuplicate(b,delta_u,ierr)
     call VecDuplicate(b,soln_prev,ierr)
     call VecDuplicate(b,temp_vec,ierr)

  !Creating vectors for local iterative solution
     !call VecCreate(PETSC_COMM_WORLD,b_local,ierr)
     !call VecSetSizes(b_local,PETSC_DECIDE,4,ierr)
     !call VecSetFromOptions(b_local,ierr)
     
     !call VecDuplicate(b_local,local_delta_u,ierr)
     !call VecDuplicate(b_local,local_soln_iter,ierr)
     !call VecDuplicate(b_local,local_soln_prev,ierr)
     !call VecDuplicate(b_local,local_b_newton,ierr)

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
  !   allocate(A_global(num_nodes,num_nodes))

     allocate(b_global(num_nodes))
     allocate(b_global_shell(num_nodes))
     allocate(qnodes(0:num_data_test%num_quadrature_nodes-1))
     allocate(qweights(0:num_data_test%num_quadrature_nodes-1))

  !   allocate(b_test(num_nodes))

     allocate(f_val(num_nodes))
  

  !   allocate(IPIV(num_nodes))

  !test_node = solve_node(1, prob_data_test, num_data)
  !test_element1 = solve_element(6, prob_data_test, num_data_test)
  !write(*,*) r(0.25_dp,test_element1)
  !write(*,*) s(0.8_dp,test_element1)
  !test_element2 = solve_element(2, prob_data_test, num_data_test)
  !test_element3 = solve_element(3, prob_data_test, num_data_test)
  !test_element4 = solve_element(4, prob_data_test, num_data_test)

  !write(*,*) test_node%neigh_elt_nos
  !write(*,*) test_element%nodes%coord(1)

     !call build_global_matrices(A_global, b_global, prob_data_test, num_data_test)

!========================Building global b===================================

     call build_global_b(b_global_shell,prob_data_test,num_data_test)

!========================Building global b===================================

  !write(*,*) MAXVAL(b_global_shell - b_global)

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

  !===========================Solution using LAPACK============================

!==============Using Petsc-ksp for solving=============
  !write(*,*) A_global

  !call VecCreateSeqWithArray(PETSC_COMM_SELF,PETSC_DECIDE,num_nodes,b_global,b,ierr)
  !call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
  !write(*,*) b_global

!=============Assigning values to vectors and Matrix===========================

  !Assigning values to RHS of linear system (Direct solution)
     call VecSetValues(b,num_nodes,col_ind,b_global_shell,INSERT_VALUES,ierr)
     call VecAssemblyBegin(b,ierr)
     call VecAssemblyEnd(b,ierr)

  !Assigning values to Matrix of linear system
  !   call MatSetValues(A,num_nodes,row_ind,num_nodes,col_ind,transpose(A_global),INSERT_VALUES,ierr)
  !   call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
  !   call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

  !Assigning values to Parameters (Iterative solution)
     call VecSet(soln_init,0.0_dp,ierr)
     call VecAssemblyBegin(soln_init,ierr)
     call VecAssemblyEnd(soln_init,ierr)
  !call VecView(soln_init,PETSC_VIEWER_STDOUT_WORLD,ierr)
  !call MatView(A,PETSC_VIEWER_STDOUT_SELF,ierr)

     call VecCopy(soln_init,soln_iter,ierr)

     call VecSet(delta_u,1.0_dp,ierr)

  !Assigning values to Parameters (Local Iterative Solution)
  !   call VecSet(local_delta_u,1.0_dp,ierr)
  
!===========Solving linear system directly=======================================

  !   call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  !   call KSPSetFromOptions(ksp,ierr)
  
  !   call KSPSetOperators(ksp,A,A,ierr)

  !   call KSPSetType(ksp,KSPCG,ierr)

  !   call KSPGetPC(ksp,pc,ierr)
  !   call PCSetType(pc,PCJACOBI,ierr)

  !do
  !     call KSPSolve(ksp,b,soln,ierr)
  !end do

  !call KSPSetFromOptions(ksp,ierr)
  

  !call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
  !write(*,*) b_global

!===========Solving linear system directly=======================================

!===========Solving linear system iteratively====================================

!     call KSPCreate(PETSC_COMM_WORLD,ksp_iter,ierr)
     call KSPCreate(PETSC_COMM_WORLD,ksp_iter_shell,ierr)
!     call KSPCreate(PETSC_COMM_WORLD,ksp_local_shell,ierr)

     call KSPSetOptionsPrefix(ksp_iter_shell,"shell_",ierr)
!     call KSPSetOptionsPrefix(ksp_local_shell,"shell_",ierr)

!     call KSPSetOperators(ksp_iter,A,A,ierr)
     call KSPSetOperators(ksp_iter_shell,A_global_shell,A_global_shell,ierr)
!     call KSPSetOperators(ksp_local_shell,A_local_shell,A_local_shell,ierr)

!     call KSPSetType(ksp_iter,KSPCG,ierr)
     call KSPSetType(ksp_iter_shell,KSPCG,ierr)
!     call KSPSetType(ksp_local_shell,KSPCG,ierr)

!     call KSPGetPC(ksp_iter,pc,ierr)
     call KSPGetPC(ksp_iter_shell,pc_shell,ierr)
!     call KSPGetPC(ksp_local_shell,pc_shell_local,ierr)

!     call PCSetType(pc,PCJACOBI,ierr)
     !call PCSetType(pc_shell,PCJACOBI,ierr)
     call PCSetType(pc_shell,PCSHELL,ierr)
!     call PCSetType(pc_shell_local,PCSHELL,ierr)

     call PCShellSetContext(pc_shell,ctxA,ierr)
!     call PCShellSetContext(pc_shell_local,ctxA,ierr)
     
!     call KSPSetFromOptions(ksp_iter,ierr)
     call KSPSetFromOptions(ksp_iter_shell,ierr)
!     call KSPSetFromOptions(ksp_local_shell,ierr)

!Testing the preconditioner operator definition

     !call VecCreate(PETSC_COMM_WORLD,temp_op_arg,ierr)
     !call VecSetSizes(temp_op_arg,PETSC_DECIDE,num_nodes,ierr)
     !call VecSetFromOptions(temp_op_arg,ierr)

     !call VecDuplicate(temp_op_arg,temp_ret_val,ierr)

     !call PC_Shell_Jacobi(pc_shell,delta_u,temp_ret_val,ierr)

     !call VecView(temp_ret_val,PETSC_VIEWER_STDOUT_WORLD,ierr)

     call PCShellSetApply(pc_shell,PC_Shell_Jacobi,ierr)
     !call PCShellSetApply(pc_shell_local,PC_Shell_Jacobi_local,ierr)

     call VecNorm(delta_u,NORM_INFINITY,norm_delta_u,ierr)
     !call VecNorm(local_delta_u,NORM_INFINITY,norm_local_delta_u,ierr)

     !Performing scaling outside the loop so that it doesn't repeat inside
     !call VecScale(b,-1.0_dp,ierr)

     tol = 1e-4

     n_iter_newton = 0

!Solving iteratively by Newton's method using assembled Global matrix
!     do while (norm_delta_u > tol)
        !Copying soln of previous iteration to soln_prev
!        call VecCopy(soln_iter,soln_prev,ierr)

        !Performing the operation b_newton = - b + A*soln_prev
!        call MatMultAdd(A,soln_prev,b,b_newton,ierr)
        !call MatMult(A,soln_prev,b_newton,ierr)

        !call VecView(b_newton,PETSC_VIEWER_STDOUT_WORLD,ierr)

        !Performing the solve for delta_u
!        call KSPSolve(ksp_iter,b_newton,delta_u,ierr)

        !call VecView(delta_u,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)

        !Calculating soln = soln_prev - delta_u
!        call VecWAXPY(soln_iter,-1.0_dp,delta_u,soln_prev,ierr)

        !Calculating norm of delta_u
!        call VecNorm(delta_u,NORM_INFINITY,norm_delta_u,ierr)

!        n_iter_newton = n_iter_newton + 1
        !write(*,*) 'Iteration number',n_iter_newton,'Error =',norm_delta_u
!     end do

     !call VecView(soln_iter,PETSC_VIEWER_STDOUT_WORLD,ierr)

     do i=1,num_nodes
        f_val(i) = 0.0_dp
        ny = num_divs_y+1-mod(i-1,num_divs_y+1)
        nx = 1 + ((i-1)/(num_divs_y+1))
        xi = hx*(nx-1)
        yi = hy*(ny-1)
        f_val(i) = f(xi,yi)
     end do

!     call VecGetValues(soln_iter,num_nodes,col_ind,b_global,ierr)

!     write(*,*) log(dble(num_divs_x)),log(maxval(abs(b_global - f_val))) 

     call VecScale(b,-1.0_dp,ierr)

     call VecSet(delta_u,1.0_dp,ierr)
     call VecNorm(delta_u,NORM_INFINITY,norm_delta_u,ierr)
     call VecCopy(soln_init,soln_iter,ierr)

     n_iter_newton = 0

!=====================Solving strategy to be used now==========================

!Solving iteratively by Newton's method using Global Shell matrix (to be used with libCEED)
     do while (norm_delta_u > tol)
        !Copying soln of previous iteration to soln_prev
        call VecCopy(soln_iter,soln_prev,ierr)

        !Performing the operation b_newton = - b + A*soln_prev
        !call MatMultAdd(A,soln_prev,b,b_newton,ierr)
        call MyMult(A_global_shell,soln_prev,b_newton,ierr)
        call VecAXPY(b_newton,1.0_dp,b,ierr)

        !call PCGetOperators(pc_shell,pshellmat,PETSC_NULL_MAT,ierr)

        !call MatGetLocalSize(pshellmat,p,q,ierr)

        !call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)

        !Performing the solve for delta_u
        call KSPSolve(ksp_iter_shell,b_newton,delta_u,ierr)

        !call VecView(b_newton,PETSC_VIEWER_STDOUT_WORLD,ierr)

        !Calculating soln = soln_prev - delta_u
        call VecWAXPY(soln_iter,-1.0_dp,delta_u,soln_prev,ierr)

        !Calculating norm of delta_u
        call VecNorm(delta_u,NORM_INFINITY,norm_delta_u,ierr)

        !write(*,*) norm_delta_u

        n_iter_newton = n_iter_newton + 1
        !write(*,*) 'Iteration number',n_iter_newton,'Error =',norm_delta_u
     end do

     call VecGetValues(soln_iter,num_nodes,col_ind,b_global,ierr)
     write(*,*) log(dble(num_divs_x)),log(maxval(abs(b_global - f_val)))

     !call VecView(soln_iter,PETSC_VIEWER_STDOUT_WORLD,ierr)

!=====================Solving strategy to be used now==========================


!Solving iteratively by Newton's method using Local Shell matrices
!     n_iter_newton = 0
!     call VecSet(local_delta_u,1.0_dp,ierr)
!     call VecSet(soln_iter,0.0_dp,ierr)

   !Loop on element numbers
!     do i=1,num_elements

   !     i=1

!        elt = solve_element(i,prob_data_test,num_data_test)
!        elt_glonos = elt%nodes%global_num

!        elt_glonos = elt_glonos - 1

        !write(*,*) elt_glonos

!        call build_local_b_vec(prob_data_test, num_data_test, elt, b_local)

!        call VecGetValues(soln_init,4,elt_glonos,local_soln,ierr)
!        call VecSetValues(local_soln_iter,4,loc_col_ind,local_soln,INSERT_VALUES,ierr)
        !call VecSet(local_soln_iter,0.0_dp,ierr)

!        call VecScale(b_local,-1.0_dp,ierr)

!        do while (norm_local_delta_u > tol)
           !Copying soln of previous iteration to soln_prev
!           call VecCopy(local_soln_iter,local_soln_prev,ierr)

           !Performing the operation b_newton = - b + A*soln_prev
!           call MyMult_local(A_local_shell,local_soln_prev,local_b_newton,ierr)
!           call VecAXPY(local_b_newton,1.0_dp,b_local,ierr)

           !Performing the solve for delta_u
!           call KSPSolve(ksp_local_shell,local_b_newton,local_delta_u,ierr)

           !Calculating soln = soln_prev - delta_u
!           call VecWAXPY(local_soln_iter,-1.0_dp,local_delta_u,local_soln_prev,ierr)

           !Calculating norm of delta_u
!           call VecNorm(local_delta_u,NORM_INFINITY,norm_local_delta_u,ierr)

!           n_iter_newton = n_iter_newton + 1

           !write(*,*) 'Iteration number',n_iter_newton,'Error =',norm_local_delta_u
           !norm_local_delta_u = 0.0_dp
!        end do

!        call VecGetValues(local_soln_iter,4,loc_col_ind,local_soln,ierr)

        !write(*,*) local_soln

!        call VecSetValues(soln_iter,4,elt_glonos,local_soln,INSERT_VALUES,ierr)
        
!        norm_local_delta_u = 1.0_dp

!     end do

     !call VecView(soln_iter,PETSC_VIEWER_STDOUT_WORLD,ierr)

!     call VecGetValues(soln_iter,num_nodes,col_ind,b_global,ierr)
!     write(*,*) log(dble(num_divs_x)),log(maxval(abs(b_global - f_val))) 

     !call VecView(soln_iter,PETSC_VIEWER_STDOUT_WORLD,ierr)


     !call VecCreate(PETSC_COMM_WORLD,temp_vec,ierr)
     !call VecSetSizes(temp_vec,PETSC_DECIDE,n,ierr)
     !call VecSetFromOptions(temp_vec,ierr)
     !call VecWAXPY(temp_vec,-1.0_dp,soln,soln_iter,ierr)
     !call VecNorm(temp_vec,NORM_INFINITY,temp_norm,ierr)
     !write(*,*) 'Infinity norm of difference between direct and iterative solution is',temp_norm

     

!===========Solving linear system iteratively====================================

!======================Plotting using Petsc Viewer===============================

  !call PetscViewerDrawGetDraw(viewer,0,draw,ierr)

!======================Plotting using Petsc Viewer===============================

     
     

!==============Deallocating Array Memory========================================

     deallocate(row_ind)
     deallocate(col_ind)
     !deallocate(A_global)
     deallocate(b_global)
     deallocate(b_global_shell)
     deallocate(qnodes)
     deallocate(qweights)
     !deallocate(b_test)
     deallocate(f_val)
     !deallocate(IPIV)

!==============Deallocating Array Memory========================================

!==============Destroying and Finalizing Petsc objects==

     !call KSPDestroy(ksp,ierr)
     !call KSPDestroy(ksp_iter,ierr)
     call KSPDestroy(ksp_iter_shell,ierr)
     !call KSPDestroy(ksp_local_shell,ierr)
     !call MatDestroy(A,ierr)
     call MatDestroy(A_global_shell,ierr)
     call VecDestroy(b,ierr)
     !call VecDestroy(soln,ierr)
     call VecDestroy(soln_iter,ierr)
     call VecDestroy(b_newton,ierr)
     call VecDestroy(soln_init,ierr)
     call VecDestroy(delta_u,ierr)
     call VecDestroy(soln_prev,ierr)
     call VecDestroy(temp_vec,ierr)
     call MatDestroy(A_local_shell,ierr)

     big_loop_variable = big_loop_variable*2

!  end do

  call PetscFinalize(ierr)

!==============Destroying and Finalizing Petsc objects==
end program main
