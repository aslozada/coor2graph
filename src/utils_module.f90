module utils_module
  use iso_fortran_env, only: output => output_unit
  implicit none

  character(len=13), parameter :: hash="(  100('-') )"

contains

  ! status = old(0)-new(1)-nuknown(2)     

  subroutine openFile(unit,buffer,status)
    integer, intent(out)     :: unit
    character(*), intent(in), optional :: buffer 
    integer, intent(in)      :: status
    ! local variables
    logical :: is
    character(len=:), allocatable :: fileName


    unit = 10
    is = .false.
    inquire(unit,opened=is)
    do 
       if(.not.is) exit
       unit = unit + 1
       inquire(unit,opened=is)
    end do

      if(status == 0) then
        if(present(buffer)) then
          fileName = trim(adjustl(buffer))

          inquire(file=fileName,exist=is)
          if(.not.is) then
            write(output,'(a)') 'No such file: '//fileName
            stop
          else
            write(output,'(a)') 'Using file: '//fileName     
          end if
        end if
      end if

    select case(status)
      case(0)
         if(present(buffer)) then
         open(unit,file=fileName,status='old')
         end if
      case(1) 
         if(present(buffer)) then
         open(unit,file=fileName,status='new')
         end if
      case(2)
         if(present(buffer)) then
         open(unit,file=fileName,status='unknown')
         end if
      case(3) 
       open(unit,status='scratch')
    end select

  end subroutine openFile

          
end module utils_module
