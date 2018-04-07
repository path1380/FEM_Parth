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
  real(kind=dp), dimension(4) :: b_local,b_local1,b_local2,b_local3,b_local4
  real(kind=dp) :: test_x,test_y

  real(kind=dp), dimension(9,9) :: A_global
  real(kind=dp), dimension(9) :: b_global

  type(problem_data) :: prob_data
  type(numerics_data) :: num_data
  type(node) :: test_node
  type(element) :: test_element,test_element1,test_element2,test_element3,test_element4

  integer :: i

!LAPACK variables

  integer, dimension(9) :: IPIV
  integer :: INFO

  !Giving problem data input
  prob_data%domain_length = 1.0_dp
  prob_data%domain_width = 1.0_dp

  !Giving numerics data input
  num_data%num_divs_x = 2
  num_data%num_divs_y = 2
  num_data%num_quadrature_nodes = 3

  !test_node = solve_node(1, prob_data, num_data)
  test_element1 = solve_element(1, prob_data, num_data)
  test_element2 = solve_element(2, prob_data, num_data)
  test_element3 = solve_element(3, prob_data, num_data)
  test_element4 = solve_element(4, prob_data, num_data)

  !write(*,*) test_node%neigh_elt_nos
  !write(*,*) test_element%nodes%coord(1)

  call build_local_A(prob_data,num_data, A_local)

  call build_global_matrices(A_global, b_global, prob_data, num_data)

  !write(*,*) A_global - transpose(A_global)
  !write(*,*) b_global
  !write(*,*) A_local - transpose(A_local)
  !write(*,*) b_local
  
  call DGETRF(9,9,A_global,9,IPIV,INFO)
  call DGETRS('N',9,1,A_global,9,IPIV,b_global,9,INFO)

  write(*,*) b_global
  
end program main
