module local_mat
#include <petsc/finclude/petsc.h>
  use type_defs
  use basic_structures
  use numbering_convention_defn
  use legendre_module
  use basis_function
  use Input_function
  use petsc
  use solver_context_interfaces
  implicit none
contains

  real(dp) function Jacobian(prob_data, num_data)
     type(numerics_data), intent(in) :: num_data
     type(problem_data), intent(in) :: prob_data
     
     real(kind=dp) :: hx,hy

     hx = (prob_data%domain_length)/(num_data%num_divs_x)
     hy = (prob_data%domain_width)/(num_data%num_divs_y)

     Jacobian = hx*hy/4.0_dp
  end function Jacobian

  subroutine build_local_A(prob_data, num_data, A_local)

!In this case, the element does not matter. In general, the Jacobian and the mapping depend on the element.
     type(problem_data) :: prob_data
     type(numerics_data) :: num_data
     real(kind=dp), dimension(4,4) :: A_local

     integer :: nq,i,j,k,l
     real(kind=dp), dimension(0:num_data%num_quadrature_nodes) :: qnodes,qweights

     nq = num_data%num_quadrature_nodes

     qnodes = num_data%nodes
     qweights = num_data%weights

     do j=1,4
        do i=1,4
           A_local(i,j) = 0.0_dp
           do k=0,nq
              do l=0,nq
                 A_local(i,j) = A_local(i,j) + (qweights(k)*qweights(l)*bilinear_basis(qnodes(k),&
                                qnodes(l),i)*bilinear_basis(qnodes(k),qnodes(l),j)*&
                                Jacobian(prob_data,num_data))
              end do
           end do
        end do
     end do
  end subroutine build_local_A


  real(kind=dp) function x(r,elt)
     real(kind=dp), intent(in) :: r
     type(element), intent(in) :: elt

     real(kind=dp) :: x1,x2
     x1 = elt%nodes(1)%coord(1)
     x2 = elt%nodes(4)%coord(1)

     x = ((x1*(1-r))+(x2*(1+r)))/2
  end function x

  real(kind=dp) function y(s,elt)
     real(kind=dp), intent(in) :: s
     type(element), intent(in) :: elt

     real(kind=dp) :: y1,y2
     y1 = elt%nodes(1)%coord(2)
     y2 = elt%nodes(2)%coord(2)

     y = ((y1*(1+s))+(y2*(1-s)))/2
  end function y

  subroutine build_local_b(prob_data, num_data, elt, b_local)
     type(problem_data) :: prob_data
     type(numerics_data) :: num_data
     type(element) :: elt
     real(kind=dp), dimension(4) :: b_local

     integer :: nq,i,j,k

     real(kind=dp), dimension(0:num_data%num_quadrature_nodes) :: qnodes,qweights

     nq = num_data%num_quadrature_nodes

     call lglnodes(qnodes,qweights,nq)

     do i=1,4
        b_local(i) = 0.0_dp
        do j=0,nq
           do k=0,nq
              b_local(i) = b_local(i) + (qweights(j)*qweights(k)*&
                           bilinear_basis(qnodes(j), qnodes(k),i)*&
                           f(x(qnodes(j),elt),y(qnodes(k),elt))*&
                           Jacobian(prob_data,num_data))
           end do
        end do
     end do
  end subroutine build_local_b

!========================Same routines returning petsc objects=====================

  subroutine build_local_A_petsc(prob_data, num_data, A_local_petsc)
     type(problem_data) :: prob_data
     type(numerics_data) :: num_data

     !Petsc Variables
     Mat A_local_petsc
     PetscScalar ierr

     integer :: nq,i,j,k,l
     real(kind=dp), dimension(0:num_data%num_quadrature_nodes) :: qnodes,qweights
     real(kind=dp) :: temp_val
 

     nq = num_data%num_quadrature_nodes

     call lglnodes(qnodes,qweights,nq)

     do j=0,3
        do i=0,3
           temp_val = 0.0_dp
           do k=0,nq
              do l=0,nq
                 temp_val = temp_val + (qweights(k)*qweights(l)*bilinear_basis(qnodes(k),&
                            qnodes(l),i+1)*bilinear_basis(qnodes(k),qnodes(l),j+1)*&
                            Jacobian(prob_data,num_data))
              end do
           end do
           call MatSetValue(A_local_petsc,i,j,temp_val,INSERT_VALUES,ierr)
        end do
     end do

  end subroutine build_local_A_petsc

  subroutine MyMult_local(local_shell, op_arg, ret_val, ierr)
     Mat local_shell
     Vec op_arg, ret_val    
     PetscErrorCode ierr

     type(MatCtx), POINTER :: local_ctx_pt

     type(problem_data) :: prob_data
     type(numerics_data) :: num_data

     integer :: nq,i,j,k,l
     integer, dimension(1) :: loc_j_array

     real(kind=dp) :: temp_val, alocalij
     real(kind=dp), dimension(1) :: temp_val_array
     real(kind=dp), allocatable, dimension(:) :: qnodes,qweights

     call MatShellGetContext(local_shell,local_ctx_pt,ierr)

     prob_data = local_ctx_pt%prob
     num_data = local_ctx_pt%num

     nq = num_data%num_quadrature_nodes

     allocate(qnodes(0:nq))
     allocate(qweights(0:nq))

     call lglnodes(qnodes,qweights,nq)

     !Initialize ret_val to 0 for adding values
     call VecSet(ret_val,0.0_dp,ierr)
     temp_val = 0.0_dp
     alocalij = 0.0_dp

     do i=1,4
        do j=1,4
           loc_j_array(1) = j-1
           do k=0,nq
              do l=0,nq
                 alocalij = alocalij + (qweights(k)*qweights(l)*bilinear_basis(qnodes(k),&
                            qnodes(l),i)*bilinear_basis(qnodes(k),qnodes(l),j)*&
                            Jacobian(prob_data,num_data))
              end do
           end do
           call VecGetValues(op_arg,1,loc_j_array,temp_val_array,ierr)
           temp_val = alocalij*temp_val_array(1)
           call VecSetValue(ret_val,i-1,temp_val,ADD_VALUES,ierr)
           alocalij = 0.0_dp
        end do
        temp_val = 0.0_dp
     end do     

     deallocate(qnodes)
     deallocate(qweights)

  end subroutine MyMult_local

  subroutine libCEED_Mult(prob_data,num_data,op_arg,ret_val)

     real(kind=dp), dimension(4) :: op_arg,ret_val
     type(problem_data) :: prob_data
     type(numerics_data) :: num_data

     real(kind=dp), dimension(4,4) :: A_local
     integer :: i,j

     call build_local_A(prob_data,num_data,A_local)

     !Initializing ret_val
     do j=1,4
        ret_val(j) = 0.0_dp
     end do

     do i=1,4
        do j=1,4
           ret_val(i) = ret_val(i) + (A_local(i,j)*op_arg(j))
        end do
     end do

  end subroutine libCEED_Mult

  subroutine build_local_b_vec(prob_data, num_data, elt, b_local_vec)
     type(problem_data) :: prob_data
     type(numerics_data) :: num_data
     type(element) :: elt
     
     Vec b_local_vec
     PetscErrorCode ierr

     integer :: nq,i,j,k

     real(kind=dp), dimension(0:num_data%num_quadrature_nodes) :: qnodes,qweights
     real(kind=dp) :: temp_val

     nq = num_data%num_quadrature_nodes

     call lglnodes(qnodes,qweights,nq)

     !Initialize the vector
     call VecSet(b_local_vec,0.0_dp,ierr)
     temp_val = 0.0_dp

     do i=1,4
        temp_val = 0.0_dp
        do j=0,nq
           do k=0,nq
              temp_val = temp_val + (qweights(j)*qweights(k)*&
                           bilinear_basis(qnodes(j), qnodes(k),i)*&
                           f(x(qnodes(j),elt),y(qnodes(k),elt))*&
                           Jacobian(prob_data,num_data))
           end do
        end do
        call VecSetValue(b_local_vec,i-1,temp_val,ADD_VALUES,ierr)
     end do
  end subroutine build_local_b_vec

  subroutine PC_Shell_Jacobi_local(pc_jacobi_local,op_arg,ret_val,ierr)

     PC pc_jacobi_local
     Vec op_arg,ret_val
     PetscErrorCode ierr

     type(MatCtx), POINTER :: local_ctx_pt

     type(problem_data) :: prob_data
     type(numerics_data) :: num_data

     integer :: nq,i,j,k,l
     integer, dimension(1) :: loc_j_array

     real(kind=dp) :: temp_val, alocalij
     real(kind=dp), dimension(1) :: temp_val_array
     real(kind=dp), allocatable, dimension(:) :: qnodes,qweights

     call PCShellGetContext(pc_jacobi_local,local_ctx_pt,ierr)

     prob_data = local_ctx_pt%prob
     num_data = local_ctx_pt%num

     nq = num_data%num_quadrature_nodes

     allocate(qnodes(0:nq))
     allocate(qweights(0:nq))

     call lglnodes(qnodes,qweights,nq)

     !Initialize ret_val to 0 for adding values
     call VecSet(ret_val,0.0_dp,ierr)
     temp_val = 0.0_dp
     alocalij = 0.0_dp

     do i=1,4
        do j=1,4
           loc_j_array(1) = j-1
           if (i == j) then
              do k=0,nq
                 do l=0,nq
                    alocalij = alocalij + (qweights(k)*qweights(l)*bilinear_basis(qnodes(k),&
                               qnodes(l),i)*bilinear_basis(qnodes(k),qnodes(l),j)*&
                               Jacobian(prob_data,num_data))
                 end do
              end do
              call VecGetValues(op_arg,1,loc_j_array,temp_val_array,ierr)
              temp_val = alocalij*temp_val_array(1)
           else if (i /= j) then
              temp_val = 0.0_dp
           end if
           call VecSetValue(ret_val,i-1,temp_val,ADD_VALUES,ierr)
           alocalij = 0.0_dp
        end do
        temp_val = 0.0_dp
     end do     

     deallocate(qnodes)
     deallocate(qweights)

  end subroutine PC_Shell_Jacobi_local

end module local_mat
