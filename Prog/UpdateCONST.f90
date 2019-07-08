!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Elementary program that read old-constants-file (possibly for older version of artemide)
!!		and create a new (for current version of artemide)
!!
!!						A.Vladimirov (08.07.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_setup
implicit none

CHARACTER(*),parameter::path=' .....    '
CHARACTER(*),parameter::fileIN='const-test'
CHARACTER(*),parameter::fileOUT='const-test+'

 call artemide_Setup_fromFile(fileIN,prefix=path)
 call CreateConstantsFile(fileOUT,prefix=path)
end program example