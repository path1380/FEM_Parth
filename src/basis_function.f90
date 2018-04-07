module basis_function
  use type_defs
  use basic_structures
  implicit none
contains
  real(kind=dp) function bilinear_basis(r,s,local_node_num)
! This function is defined in (r,s) space, from which the element is mapped to the (x,y) space. The mapping is used
! to simplify integrals over the element and all basis functions have the same form in the pre-image of the element
! in the (r,s) space. Only the Jacobian can differ from element to element, which, in this case, is not happening
! either.
     real(dp), intent(in) :: r,s
     integer, intent(in) :: local_node_num

     integer :: int1,int2

     int1 = (local_node_num-1)/2
     int2 = mod(local_node_num-1,2)

     bilinear_basis = (1 + (((-1)**int2)*s))*(1 - (((-1)**int1)*r))/4.0_dp

  end function bilinear_basis
end module basis_function
