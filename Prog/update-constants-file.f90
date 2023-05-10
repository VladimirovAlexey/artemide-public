!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	The program that check the constants-file and update it if nessecary
!!		requares an argument (full path to constants file)
!!
!!						A.Vladimirov (05.11.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_setup
use IO_functions
implicit none
  real*8::Q
  CHARACTER(len=516)::path,path_backup
  CHARACTER(100)::line
  logical::file_exists,save_backup,to_update,backup_replace
  CHARACTER(len=1)::ask
  integer::iostatus
 
 write(*,*) "--------------------------------------------------------------------------------"
 write(*,*) "-------------------- artemide.UpdateConstantsFile ------------------------------"
 
 to_update=.false.
 save_backup=.true.
 backup_replace=.false.
 
 !!! receive the argument
 call get_command_argument(1,path)
 if(len_trim(path)==0) then
  write(*,*) color(" artemide.UpdateConstantsFile-routine: pass the path to constants-file.",c_red)
  call exit(0)
 end if
 
 !!chech the presence of the file
 INQUIRE(FILE=trim(path), EXIST=file_exists) 
 if(.not.file_exists) then
  write(*,*) color(" artemide.UpdateConstantsFile-routine: constants-file is not found.",c_red)
  write(*,*) color(" check path: "//trim(path),c_red)
  call exit(0)
 end if
 
 !! output with the path
 write(*,*) " FILE: ",color(trim(path),c_green_bold)
 
 !! check for nessecity to update
 if(CheckConstantsFile(trim(path))) then
  write(*,*) color("  ... update is not required",c_green)
  call exit(0)
 else
  write(*,*) color("  ... update is required",c_green)
 end if
 
 !!ask user
 write(*,*) color("  Do you want to update this file?  (y/n)",c_green_bold)
 read(*,*) ask
 if(ask=='y') to_update=.true.
 if (.not.to_update) call exit(0)
 
 write(*,*) color("  Do you want to save backup copy of this file?  (y/n)",c_green_bold)
 read(*,*) ask
 if(ask=='n') save_backup=.false.
 
 !!! check for possibility to backup
 if(save_backup) then
  path_backup=trim(path)//"_backup"
  write(*,*) '  backup copy will be saved to :',color(trim(path_backup),c_green_bold)
  INQUIRE(FILE=trim(path_backup), EXIST=file_exists) 
    if(file_exists) then
      write(*,*) color('   backup file with such name already exists. Replace? (y/n)',c_red_bold)
      read(*,*) ask
      if(ask=='y') backup_replace=.true.
      if(.not.backup_replace) then 
	  write(*,*) '   program terminated'
	  call exit(0)
      end if
      if(backup_replace) write(*,*) '   ... backup will be replaced'
    end if
  end if
  
  !!! save backup
  if(save_backup) then
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    OPEN(UNIT=52, FILE=path_backup, ACTION="write", STATUS="replace")
    do
      read(51,'(A)',IOSTAT=iostatus) line
      if(iostatus>0) then
	write(*,*) color('    ... some problem in reading file...',c_red_bold)
	stop
      else if(iostatus<0) then
	exit
      else
	write(52,'(A)') trim(line)
      end if
    end do
    CLOSE (51, STATUS='KEEP')
    CLOSE (52, STATUS='KEEP')
    
    write(*,*) color('   ... backup file saved normally',c_green)
  end if
 
  call artemide_Setup_fromFile(trim(path))
  call CreateConstantsFile(trim(path))
  
  write(*,*) color('  DONE!',c_green) 
end program example
