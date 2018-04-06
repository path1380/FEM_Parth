!================================================================================================================
!
! File: main.f90
! Brief description: Main program of the Project.
!
! Detailed description: This program will be regulated by a perl script which gives inputs to a template file. 
!                       The code is documented in a point-wise form. Please refer to the points in order to 
!                       understand what the heck is going on in the program.
!
! PLEASE MAINTAIN THE DOCUMENTATION FORMAT OF THE CODE. IT IS ABSOLUTELY CRITICAL FOR OTHERS TO UNDERSTAND WHAT
! YOU'VE DONE.
!
! Authors: Parth Thakkar, Fortino Garcia, Alex Rybchuk, William Boshell
!
!================================================================================================================

program main
  use type_defs
  implicit none


  write(*,*) -4*ATAN(-1.0_dp)
  write(*,*) exp(1.0_dp)
end program main
