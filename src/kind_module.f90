!-------------------------------------------------------------------
! This module is part of coor2graph
! Declare the numeric kinds parameters
!-------------------------------------------------------------------
!> author: Asdrubal Lozada
!  update: 4/11/2026
!  license: GPL-v3 
  
module kind_module 

use iso_fortran_env, only: real128,real64,real32,int32
implicit none

private 

integer,parameter,public :: qp = real128
integer,parameter,public :: wp = real64
integer,parameter,public :: sp = real32
integer,parameter,public :: ip = int32

end module kind_module
