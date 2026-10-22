! SPDX-License-Identifier: MIT
!
! This file is part of MONKI (Monte Carlo KNMI), version 1.0.1.
!
! Author, creator, and lead developer:
!   Victor J.H. Trees
!
! Contact:
!   victor.trees@knmi.nl
!   v.j.h.trees@tudelft.nl
!
! Copyright (c) 2026 Victor J.H. Trees

  
MODULE write_messages
  !
  INTEGER,PRIVATE,PARAMETER    :: unit=6
  CHARACTER(len=12),PRIVATE,PARAMETER   :: debug_str='Debug    | '
  CHARACTER(len=12),PRIVATE,PARAMETER   :: error_str='Error    | '
  CHARACTER(len=12),PRIVATE,PARAMETER   ::  info_str='Info     | '
  CHARACTER(len=12),PRIVATE,PARAMETER   ::  prog_str='Progress | '
  CHARACTER(len=12),PRIVATE,PARAMETER   ::  warn_str='Warning  | '
  !
CONTAINS
  !
  SUBROUTINE write_debug(message)
    !    
    CHARACTER(len=*),INTENT(in)   :: message
    WRITE(unit,'(a)') debug_str//TRIM(ADJUSTL(message))
  END SUBROUTINE write_debug
  !
  SUBROUTINE write_error(message)
    !    
    CHARACTER(len=*),INTENT(in)   :: message
    WRITE(unit,'(a)') error_str//TRIM(ADJUSTL(message))
  END SUBROUTINE write_error
  !
  SUBROUTINE write_info(message)
    !    
    CHARACTER(len=*),INTENT(in)   :: message
    WRITE(unit,'(a)') info_str//TRIM(ADJUSTL(message))
  END SUBROUTINE write_info
  !
  SUBROUTINE write_prog(i,itot)
    !    
    Integer,intent(in)   :: i,itot
    !
    Character(len=30)    :: message
    !
    write(message,'(i6,a,i6)') i,' of ',itot
    !
    WRITE(unit,'(a)') prog_str//TRIM(ADJUSTL(message))
  END SUBROUTINE write_prog
  !
  SUBROUTINE write_prog_na(i,itot)
    !    
    Integer,intent(in)   :: i,itot
    !
    Character(len=30)    :: message
    !
    write(message,'(i6,a,i6,a1)') i,' of ',itot,achar(13)
    !
    WRITE(unit,'(a)',advance="NO") prog_str//TRIM(ADJUSTL(message))
    !
    if (i==itot) write(unit,*) ' '
    
  END SUBROUTINE write_prog_na
  !
  !
  SUBROUTINE write_warning(message)
    !    
    CHARACTER(len=*),INTENT(in)   :: message
    WRITE(unit,'(a)') warn_str//TRIM(ADJUSTL(message))
  END SUBROUTINE write_warning
  !
END MODULE write_messages
