module global_mat
#include <petsc/finclude/petsc.h>
  use type_defs
  use basic_structures
  use numbering_convention_defn
  use legendre_module
  use basis_function
  use Input_function
  use petsc
  use local_mat
  use solver_context_interfaces
  implicit none
contains
  subroutine build_global_matrices(A_global, b_global, prob_data, num_data)
     type(problem_data) :: prob_data
     type(numerics_data) :: num_data
     real(kind=dp), dimension(((num_data%num_divs_x)+1)*((num_data%num_divs_y)+1),&
                             ((num_data%num_divs_x)+1)*((num_data%num_divs_y)+1)) :: A_global
     real(kind=dp), dimension(((num_data%num_divs_x)+1)*((num_data%num_divs_y)+1)) :: b_global

     integer :: num_divs_x,num_divs_y,num_elements,num_nodes,num_quadrature_nodes
     integer :: i,j,k,l,glo_i,glo_j

     type(element) :: temp_element
     real(kind=dp), dimension(4,4) :: temp_A_local
     real(kind=dp), dimension(4) :: temp_b_local

     num_divs_x = num_data%num_divs_x
     num_divs_y = num_data%num_divs_y
     num_quadrature_nodes = num_data%num_quadrature_nodes
     
     num_elements = num_divs_x*num_divs_y
     num_nodes = (num_divs_x+1)*(num_divs_y+1)

     !Initialize the global matrices

     do j=1,num_nodes
        do i=1,num_nodes
           A_global(i,j) = 0.0_dp
        end do
        b_global(j) = 0.0_dp
     end do

     !Build the global matrices
     do k=1,num_elements
        temp_element = solve_element(k,prob_data,num_data)

        call build_local_A(prob_data,num_data,temp_A_local)
        call build_local_b(prob_data,num_data,temp_element,temp_b_local)

!glo_i and glo_j are the global node numbers of nodes in element temp_element with local node numbers i and j.

        do j=1,4
           glo_j = temp_element%nodes(j)%global_num
           do i=1,4
              glo_i = temp_element%nodes(i)%global_num
              A_global(glo_i,glo_j) = A_global(glo_i,glo_j) + temp_A_local(i,j)
           end do
           b_global(glo_j) = b_global(glo_j) + temp_b_local(j)
        end do
     end do
  end subroutine build_global_matrices

  subroutine MyMult(global_shell,op_arg,ret_val,ierr)

!======================================================================================
!
!  Brief description : This subroutine defines the global matrix as an operator. 
!
!  Inputs : global_shell - Shell matrix for which operator is to be defined
!           op_arg       - Vector on which the operation is to be carried out
!           ierr         - Error code
!
!  Output : ret_val      - Vector obtained as a result of the operation
!
!  Detailed description : The subroutine takes the shell matrix (global_shell) as input 
!                         and defines it as an operator on the operator argument (op_arg)
!                         which is equivalent to matrix multiplication to return the 
!                         result (ret_val)
!
!======================================================================================
     
     Mat global_shell
     Vec op_arg,ret_val
  
     PetscErrorCode ierr

     type(MatCtx), POINTER :: ctx_pt

     type(problem_data) :: prob_data
     type(numerics_data) :: num_data

     integer :: num_divs_x,num_divs_y,num_elements,num_nodes,num_quadrature_nodes
     integer :: i,j,k,l,glo_i,glo_j
     integer, dimension(1) :: glo_j_array

     real(kind=dp), dimension(4,4) :: temp_A_local

     !Value in glo_ith row and glo_jth column of matrix represented by global_shell
     real(kind=dp) :: global_shell_gloigloj,temp_val
     real(kind=dp), dimension(1) :: temp_val_array

     type(element) :: temp_element
     
     call MatShellGetContext(global_shell,ctx_pt,ierr)
     temp_A_local = ctx_pt%local_matrix
     prob_data = ctx_pt%prob
     num_data = ctx_pt%num

     num_divs_x = num_data%num_divs_x
     num_divs_y = num_data%num_divs_y
     num_quadrature_nodes = num_data%num_quadrature_nodes
     
     num_elements = num_divs_x*num_divs_y
     num_nodes = (num_divs_x+1)*(num_divs_y+1)

     global_shell_gloigloj = 0.0_dp

     do k=1,num_elements
        temp_element = solve_element(k,prob_data,num_data)

!        call build_local_A(prob_data,num_data,temp_A_local)
!        call build_local_b(prob_data,num_data,temp_element,temp_b_local)

!glo_i and glo_j are the global node numbers of nodes in element temp_element with local node numbers i and j.

        do j=1,4
           glo_j = temp_element%nodes(j)%global_num
           glo_j_array(1) = glo_j-1
           do i=1,4
              glo_i = temp_element%nodes(i)%global_num
              global_shell_gloigloj = temp_A_local(i,j)
              call VecGetValues(op_arg,1,glo_j_array,temp_val_array,ierr)
              temp_val = global_shell_gloigloj*temp_val_array(1)
              call VecSetValue(ret_val,glo_i-1,temp_val,ADD_VALUES,ierr)
!              ret_val(glo_i) = ret_val(glo_i) + (global_shell_gloigloj*op_arg(glo_j))
           end do
!           b_global(glo_j) = b_global(glo_j) + temp_b_local(j)
        end do
     end do

  end subroutine MyMult

  subroutine build_global_b(b_global, prob_data, num_data)
     type(problem_data) :: prob_data
     type(numerics_data) :: num_data
     
     real(kind=dp), dimension(((num_data%num_divs_x)+1)*((num_data%num_divs_y)+1)) :: b_global

     integer :: num_divs_x,num_divs_y,num_elements,num_nodes,num_quadrature_nodes
     integer :: i,j,k,l,glo_i,glo_j

     type(element) :: temp_element
     
     real(kind=dp), dimension(4) :: temp_b_local

     num_divs_x = num_data%num_divs_x
     num_divs_y = num_data%num_divs_y
     num_quadrature_nodes = num_data%num_quadrature_nodes
     
     num_elements = num_divs_x*num_divs_y
     num_nodes = (num_divs_x+1)*(num_divs_y+1)

     !Initialize the global matrices

     do j=1,num_nodes
        b_global(j) = 0.0_dp
     end do

     !Build the global matrices
     do k=1,num_elements
        temp_element = solve_element(k,prob_data,num_data)

        call build_local_b(prob_data,num_data,temp_element,temp_b_local)

!glo_i and glo_j are the global node numbers of nodes in element temp_element with local node numbers i and j.

        do j=1,4
           glo_j = temp_element%nodes(j)%global_num
           b_global(glo_j) = b_global(glo_j) + temp_b_local(j)
        end do
     end do     
  end subroutine build_global_b

  subroutine PC_Shell_Jacobi(pc_jacobi,op_arg,ret_val,ierr)
!======================================================================================
!
!  Brief description : This subroutine defines the preconditioner as an operator. 
!
!  Inputs : pc_jacobi    - Shell Preconditioner for which operator is to be defined
!           op_arg       - Vector on which the operation is to be carried out
!           ierr         - Error code
!
!  Output : ret_val      - Vector obtained as a result of the operation
!
!  Detailed description : The subroutine takes the shell preconditioner (pc_jacobi)
!                         as input and defines it as an operator on the operator 
!                         argument (op_arg) which is equivalent to matrix 
!                         multiplication to return the result (ret_val)
!
!======================================================================================

     PC pc_jacobi
     Vec op_arg,ret_val
     PetscErrorCode ierr

     type(MatCtx), POINTER :: ctx

     type(problem_data) :: prob_data
     type(numerics_data) :: num_data

     integer :: num_divs_x,num_divs_y,num_elements,num_nodes,num_quadrature_nodes
     integer :: i,j,k,l,glo_i,glo_j
     integer, dimension(1) :: glo_j_array

     real(kind=dp), dimension(4,4) :: temp_A_local

     !Value in glo_ith row and glo_jth column of matrix represented by global_shell
     real(kind=dp) :: global_shell_gloigloj,temp_val
     real(kind=dp), dimension(1) :: temp_val_array

     type(element) :: temp_element
     
     call PCShellGetContext(pc_jacobi,ctx,ierr)

     temp_A_local = ctx%local_matrix
     prob_data = ctx%prob
     num_data = ctx%num

     num_divs_x = num_data%num_divs_x
     num_divs_y = num_data%num_divs_y
     num_quadrature_nodes = num_data%num_quadrature_nodes
     
     num_elements = num_divs_x*num_divs_y
     num_nodes = (num_divs_x+1)*(num_divs_y+1)

     global_shell_gloigloj = 0.0_dp

     do k=1,num_elements
        temp_element = solve_element(k,prob_data,num_data)

!        call build_local_A(prob_data,num_data,temp_A_local)
!        call build_local_b(prob_data,num_data,temp_element,temp_b_local)

!glo_i and glo_j are the global node numbers of nodes in element temp_element with local node numbers i and j.

        do j=1,4
           glo_j = temp_element%nodes(j)%global_num
           glo_j_array(1) = glo_j-1
           do i=1,4
              glo_i = temp_element%nodes(i)%global_num
              if (glo_i /= glo_j) then
                 global_shell_gloigloj = 0.0_dp
              else if (glo_i == glo_j) then
                 global_shell_gloigloj = temp_A_local(i,j)
              end if
              call VecGetValues(op_arg,1,glo_j_array,temp_val_array,ierr)
              temp_val = global_shell_gloigloj*temp_val_array(1)
              call VecSetValue(ret_val,glo_i-1,temp_val,ADD_VALUES,ierr)
!              ret_val(glo_i) = ret_val(glo_i) + (global_shell_gloigloj*op_arg(glo_j))
           end do
!           b_global(glo_j) = b_global(glo_j) + temp_b_local(j)
        end do
     end do

  end subroutine PC_Shell_Jacobi
end module global_mat
