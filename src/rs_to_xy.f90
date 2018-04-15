module rs_to_xy
   use type_defs
   use basic_structures
   use numbering_convention_defn
   implicit none
contains
   real(kind=dp) function r(x,el)
      type(element), intent(in) :: el
      real(kind=dp), intent(in) :: x
      
      real(kind=dp) :: x0,x1

      x0 = el%nodes(1)%coord(1)
      x1 = el%nodes(3)%coord(1)

      r = ((2*x)-x0-x1)/(x1-x0)
   end function r

   real(kind=dp) function s(y,el)
      type(element), intent(in) :: el
      real(kind=dp), intent(in) :: y

      real(kind=dp) :: y0,y1

      y0 = el%nodes(2)%coord(2)
      y1 = el%nodes(1)%coord(2)

      s = ((2*y)-y0-y1)/(y1-y0)
   end function s

end module rs_to_xy
