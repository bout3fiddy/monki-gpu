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
! Copyright (c) 2026 Victor J.H. Trees!
!-----------------------------------------------------------------------

SUBROUTINE get_O3_crabs(monki_home,absfile,wav,nz,temp,crabs)
!************************************************************************
!* Get the absorption cross-sections of O3 gas
!* Data is taken from O3T_BP.dat and interpolated to the input
!* wavelength, and to the temperatures in the layers. 
!*
!* Victor Trees, 2 Jan 2023
!************************************************************************

   IMPLICIT NONE

   INTEGER,INTENT(IN) :: nz
   CHARACTER*72,INTENT(IN) :: absfile
   REAL*8,INTENT(IN) :: temp(nz+1),wav
   CHARACTER(len=300),INTENT(IN) :: monki_home
   REAL*8,INTENT(OUT) :: crabs(nz+1)

   INTEGER :: nclines,nabswavs,ntempcoefs,i,j,itemp,ilev

   REAL*8 :: tempC,temp_offset

   REAL*8,DIMENSION(:),ALLOCATABLE :: abswav,crabsint

   REAL*8,DIMENSION(:,:),ALLOCATABLE :: crabstab
   
   CHARACTER*72 :: str
   CHARACTER*200 :: absfilename

!------------------------------------------------------------------------
!  Settings
!------------------------------------------------------------------------  
   ntempcoefs = 3
   temp_offset=0.D0

!------------------------------------------------------------------------
!  Open the file
!------------------------------------------------------------------------
   absfilename = TRIM(ADJUSTL(monki_home))//TRIM(ADJUSTL(absfile))
   OPEN(file=absfilename,unit=91,err=900,status='old')

!------------------------------------------------------------------------
!  Read the number of comment lines and number of wavelengths
!------------------------------------------------------------------------   
   READ(91,*,end=999) nclines
   READ(91,*,end=999) nabswavs

   !PRINT *,'nclines=',nclines
   !PRINT *,'nabswavs=',nabswavs

!------------------------------------------------------------------------
!  Read the comment lines
!------------------------------------------------------------------------   
   DO i=1,nclines
      READ(91,2,end=999) str
      !PRINT *,str
   ENDDO

!------------------------------------------------------------------------
!  Allocate and initialize gas absorption cross-section arrays
!------------------------------------------------------------------------
   ALLOCATE(crabstab(nabswavs,ntempcoefs))
   ALLOCATE(abswav(nabswavs))
   ALLOCATE(crabsint(ntempcoefs))

   DO i=1,nabswavs
      DO j=1,ntempcoefs
         crabstab(i,j)=0.d0
      ENDDO
   ENDDO

!------------------------------------------------------------------------
!  Read in the lines containing the absorption cross-sections.
!------------------------------------------------------------------------   
   DO i=1,nabswavs
      READ(91,*,end=999) abswav(i),(crabstab(i,itemp),itemp=1,ntempcoefs)
      !READ(91,*,end=999) abswav(i),crabstab(i,1),crabstab(i,2),crabstab(i,3)
      !PRINT *,abswav(i),crabstab(i,1),crabstab(i,2),crabstab(i,3)
      !PRINT *,'abswav(i)=',abswav(i)
      !PRINT *,'crabstab(i,:)=',crabstab(i,:)
   ENDDO

!------------------------------------------------------------------------
!  Close the file
!------------------------------------------------------------------------ 
   CLOSE(unit=91)

!------------------------------------------------------------------------
!  Linear interpolation in absorption cross-section table
!  to find absorption cross-section at *waveln*, *crabs*.
!------------------------------------------------------------------------     
   DO itemp=1,ntempcoefs
      CALL linpol_dp_single(abswav,crabstab(:,itemp),nabswavs,wav,crabsint(itemp))
   ENDDO

   !PRINT *,'crabsint=',crabsint

!------------------------------------------------------------------------ 
!   Evaluate the temperature-dependence relation:
!
!    crabs = crabs0 + crabs1*temperature + crabs2*(temperature**2)
!
!        *crabs* contains the abs. cross-section for the 
!        temperature for all levels.
!        The atmospheric temperature profile is in Kelvin.
!        The temperature dependence is converted to Celcius, using
!        the parameter temp_offset.
!------------------------------------------------------------------------ 
   DO ilev=1,nz+1

      crabs(ilev) = 0.D0
      tempC = temp(ilev) - 273.15D0

      DO itemp=1,ntempcoefs

         crabs(ilev) = crabs(ilev) + &
     &               crabsint(itemp)*((tempC-temp_offset)**(itemp-1))

      ENDDO

   ENDDO
      

   !PRINT *,'crabs =',crabs
   !READ(*,*)

!------------------------------------------------------------------------
!  Deallocate
!------------------------------------------------------------------------
   DEALLOCATE(crabstab)
   DEALLOCATE(abswav)
   DEALLOCATE(crabsint)

!------------------------------------------------------------------------   
   RETURN

!------------------------------------------------------------------------
2   FORMAT(a72)
900 PRINT *,'absfilename=',absfilename
    STOP 'get_O3_crabs: error opening the absorption file'
999 PRINT *,'get_O3_crabs: unexpected end of file '
    STOP 'in get_O3_crabs unexpected end of file encountered'

END SUBROUTINE get_O3_crabs
