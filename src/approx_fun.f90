module approx_fun
  use type_defs
  use basic_structures
  use basis_function
  use numbering_convention_defn
  use rs_to_xy
  implicit none 
contains
  real(kind=dp) function approx_eval(x,y,u_l2,prob_data,num_data)
     real(kind=dp), intent(in) :: x,y
     type(problem_data), intent(in) :: prob_data
     type(numerics_data), intent(in) :: num_data
     real(kind=dp), dimension((num_data%num_divs_x+1)*(num_data%num_divs_y+1)), intent(in) :: u_l2

     type(element) :: el

     real(kind=dp) :: length,width
     real(kind=dp) :: hx,hy

     integer :: num_divs_x,num_divs_y

     integer :: num_elements,el_no_x,el_no_y,el_no
     integer :: i

     length = prob_data%domain_length
     width = prob_data%domain_width

     num_divs_x = num_data%num_divs_x
     num_divs_y = num_data%num_divs_y
     

     num_elements = num_divs_x*num_divs_y
     
     hx = length/num_divs_x
     hy = width/num_divs_y

     el_no_x = 1 + int(x/hx)
     el_no_y = 1 + int(y/hy)

     if (el_no_x == num_divs_x+1) then
        el_no_x = num_divs_x
     end if

     if (el_no_y == num_divs_y+1) then
        el_no_y = num_divs_y
     end if

     el_no = (num_divs_y - el_no_y + 1) + ((el_no_x - 1)*num_divs_y)

     el = solve_element(el_no,prob_data,num_data)

     approx_eval = 0.0_dp

     do i=1,4
        approx_eval = approx_eval + (u_l2(el%nodes(i)%global_num)*bilinear_basis(r(x,el),s(y,el),i))
     end do

     write(*,*) el_no_x,el_no_y,el_no
  end function approx_eval
end module approx_fun
