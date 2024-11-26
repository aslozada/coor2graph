!----------------------------------------------------------------------------------------------------------
!> COOR2GRAPH - Software to build an adjancency matrix from a coordinate file and compute a graph measure
!
!> author: Asdrubal Lozada
!> e-mail: aslozada@gmail.com
!----------------------------------------------------------------------------------------------------------
program main
  use kind_module
  use utils_module
  use commands_module
  use frame_module
  use networkx_module

  integer :: unit
  type(frame) :: gro

  call commands()
  call openFile(unit,groFile,0)

  call on_the_fly(gro,unit)

end program main
