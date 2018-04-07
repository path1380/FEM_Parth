module Input_function
   use type_defs
   implicit none
contains
   real(kind=dp) function f(x,y)
      real(kind=dp), intent(in) :: x,y
      f = x
      !f = x*y*(1.0_dp-x)*(1.0_dp-y)
   end function f
end module Input_function
