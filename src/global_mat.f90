module global_mat
  use type_defs
  use basic_structures
  use numbering_convention_defn
  use legendre_module
  use basis_function
  use Input_function
  use local_mat
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
end module global_mat
