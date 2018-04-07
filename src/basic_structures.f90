module basic_structures
  use type_defs
  implicit none

  type node
     integer :: global_num
     integer, dimension(4) :: neigh_elt_nos
     integer, dimension(4) :: local_num_in_elt
     real(kind=dp), dimension(2) :: coord
  end type node

  type element
     integer :: element_num
     type(node), dimension(4) :: nodes
  end type element

  type problem_data
     real(kind=dp) :: domain_length,domain_width
  end type problem_data

  type numerics_data
     integer :: num_divs_x,num_divs_y,num_quadrature_nodes
  end type numerics_data
!contains
  
end module basic_structures
