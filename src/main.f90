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
  use type_defs
  use basic_structures
  use basis_function
  use Input_function
  use numbering_convention_defn
  use legendre_module
  use local_mat
  use global_mat
  implicit none

  real(kind=dp), dimension(0:2) :: qnodes,qweights
  real(kind=dp), dimension(4,4) :: A_local
  real(kind=dp), dimension(4) :: b_local
  real(kind=dp) :: test_x,test_y

  type(problem_data) :: prob_data
  type(numerics_data) :: num_data
  type(node) :: test_node
  type(element) :: test_element

  !Giving problem data input
  prob_data%domain_length = 1.0_dp
  prob_data%domain_width = 1.0_dp

  !Giving numerics data input
  num_data%num_divs_x = 6
  num_data%num_divs_y = 5
  num_data%num_quadrature_nodes = 3

  test_node = solve_node(1, prob_data, num_data)
  test_element = solve_element(6, prob_data, num_data)

  write(*,*) test_node%neigh_elt_nos
  write(*,*) test_element%nodes%coord(1)

  call build_local_A(prob_data,num_data, A_local)
  call build_local_b(prob_data,num_data, test_element, b_local)
  
end program main
