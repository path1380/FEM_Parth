module local_mat
  use type_defs
  use basic_structures
  use numbering_convention_defn
  use legendre_module
  use basis_function
  use Input_function
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

     call lglnodes(qnodes,qweights,nq)

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
end module local_mat
