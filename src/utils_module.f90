!-------------------------------------------------------------------
! This module is part of coor2graph
!-------------------------------------------------------------------
!> author: Asdrubal Lozada
!  update: 4/11/2026
!  license: GPL-v3 
  
module utils_module

use iso_fortran_env, only: error_unit
implicit none

private

character(len=13), parameter :: hash="(  100('-') )"

public :: hash, openFile, fail
contains

  ! status = 'old', 'new', 'replace', 'scratch'  
  subroutine openFile(unit,buffer,status)
    integer, intent(out)               :: unit
    character(*), intent(in), optional :: buffer 
    character(*), intent(in)           :: status
    character(len=:), allocatable      :: name
    integer :: ios

    if (status == 'scratch') then
        open(newunit=unit, status='scratch', iostat=ios)
        if (ios /=0 ) call fail('cannot a open a scratch file')
        return
    end if

    if (.not. present(buffer)) call fail('openFile: a file name is required')
    name = trim(adjustl(buffer))
    if (len(name) == 0) call fail('openFile: the file name is empty')

    select case(status)
      case('old')
         open(newunit=unit, file=name, status='old', action='read', iostat=ios)
      case('new') 
         open(newunit=unit, file=name, status='new', action='write', iostat=ios)
      case('replace')
         open(newunit=unit, file=name, status='replace', action='write', iostat=ios)
      case default 
         call fail('openFile: unknown status: '//status)
    end select

    if(ios /= 0 ) call fail('cannot open the file: '//name)

  end subroutine openFile

  subroutine fail(msg)
    character(*), intent(in) :: msg
    
    write(error_unit,'(a)') 'coor2graph: error: '//msg
    stop 1
  end subroutine fail

          
end module utils_module
