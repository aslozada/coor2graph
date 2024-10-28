!------------------------------------------------------
program main
  use kind_module
  use utils_module
  use commands_module
  use frame_module
  use networkx_module
  ! local variable
  integer :: unit

  type(frame) :: gro


  call commands()
  call openFile(unit,fileInput,0)

  ! test ising
  !@ call ising()
  !-------------------------

  call on_the_fly(gro,unit)

  ! test

  call build_script()
!  call gro%read(unit)

end program main
