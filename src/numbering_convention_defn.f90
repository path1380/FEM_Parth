module numbering_convention_defn
  use type_defs
  use basic_structures
contains
  type(node) function solve_node(global_node_num, prob_data, num_data)
! This function solves for all node parameters just taking the global node number as input!
     integer, intent(in) :: global_node_num
     type(problem_data), intent(in) :: prob_data
     type(numerics_data), intent(in) :: num_data

     integer :: node_num_x, node_num_y
     integer :: num_divs_x, num_divs_y
     real(dp) :: hx,hy
     integer, dimension(4) :: neigh_elt_nos, local_num_in_elt

     num_divs_x = num_data%num_divs_x
     num_divs_y = num_data%num_divs_y

     node_num_x = ((global_node_num-1)/(num_divs_y+1))+1
     node_num_y = num_divs_y + 1 - mod(global_node_num-1,num_divs_y+1)

     hx = (prob_data%domain_length)/num_divs_x
     hy = (prob_data%domain_width)/num_divs_y

     neigh_elt_nos(1) = ((node_num_x-2)*num_divs_y) + (num_divs_y + 1 - node_num_y)
     neigh_elt_nos(2) = neigh_elt_nos(1) + 1
     neigh_elt_nos(3) = ((node_num_x-1)*num_divs_y) + (num_divs_y + 1 - node_num_y)
     neigh_elt_nos(4) = neigh_elt_nos(3) + 1

     !Only the appropriate node numbers will be accessed in the implementation.

     local_num_in_elt = (/4,3,2,1/)
     

     solve_node%global_num = global_node_num
     solve_node%coord = (/(node_num_x-1)*hx,(node_num_y-1)*hy/)
     solve_node%neigh_elt_nos = neigh_elt_nos
     solve_node%local_num_in_elt = local_num_in_elt
  end function solve_node

  type(element) function solve_element(elt_num, prob_data, num_data)
     integer, intent(in) :: elt_num
     type(problem_data), intent(in) :: prob_data
     type(numerics_data), intent(in) :: num_data

     integer :: elt_num_x, elt_num_y
     integer :: num_divs_x, num_divs_y
     integer :: i
     integer, dimension(4) :: global_node_nums(4)
     type(node), dimension(4) :: temp_node

     num_divs_x = num_data%num_divs_x
     num_divs_y = num_data%num_divs_y

     elt_num_x = ((elt_num-1)/num_divs_y) + 1
     elt_num_y = mod(elt_num-1,num_divs_y) + 1

     
     global_node_nums(1) = ((elt_num_x-1)*(num_divs_y+1))+(elt_num_y)
     global_node_nums(2) = global_node_nums(1) + 1
     global_node_nums(3) = (elt_num_x*(num_divs_y+1))+(elt_num_y)
     global_node_nums(4) = global_node_nums(3) + 1

     do i=1,4
        temp_node(i) = solve_node(global_node_nums(i), prob_data, num_data)
     end do

     solve_element%element_num = elt_num
     solve_element%nodes = temp_node

     !write(*,*) elt_num_x, elt_num_y
  end function solve_element
end module numbering_convention_defn
