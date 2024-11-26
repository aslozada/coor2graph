!-------------------------------------------------------------------
!> author: Asdrubal Lozada
!
! Define the numeric kinds
  module kind_module 

    use iso_fortran_env, only: real64,real32,int32,output_unit

    private 

    integer,parameter,public :: wp = real64
    integer,parameter,public :: sp = real32
    integer,parameter,public :: ip = int32
    integer,parameter,public :: stdout = output_unit

  end module kind_module
