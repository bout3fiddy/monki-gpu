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
!-----------------------------------------------------------------------
! MODULE main_routines_cartesian
!
! Module that contains subroutines for the main.f90 program of MONKI
! Subroutines that apply to the Cartesian version of the code.
!
! Author: Victor Trees
!-------------------------------------------------------------------
MODULE main_routines_cartesian

  IMPLICIT NONE

CONTAINS

!+-----------------------------------------------------------------
! SUBROUTINE rdsize
!+-----------------------------------------------------------------
! Read the scene dimensions from the input file  
!+-----------------------------------------------------------------
  SUBROUTINE rdsize(infile,nx,ny,nz,x_dim,y_dim,nmat,nscatt_types,&
    &         n_scat_angles,wave_type,wave_min,wave_max,dwave,ngases,nclouds)

    IMPLICIT NONE
    
    ! Input -------------------------------------------------------
    CHARACTER(len=500),INTENT(IN) :: infile   !Name of input file

    ! Output ------------------------------------------------------
    INTEGER,INTENT(OUT)      :: nx,ny,nz      !Number of grid cells in x, y and z directions (-)
    INTEGER,INTENT(OUT)      :: nmat          !Size of the Stokes vector: 1, 3 or 4 (-)
    INTEGER,INTENT(OUT)      :: nscatt_types  !Number of cloud droplet types (-)
    INTEGER,INTENT(OUT)      :: ngases        !Number of absorbing gases (-)
    INTEGER,INTENT(OUT)      :: nclouds       !Maximum number of different clouds in a layer (-)
    INTEGER,INTENT(OUT)      :: n_scat_angles !Number of scattering angles in Mie LUTs (-)

    REAL*8,INTENT(OUT)       :: x_dim,y_dim   !Dimensions of the scene in x and y directions (km)
    REAL*8,INTENT(OUT)       :: wave_min      !Minimum wavelength or wavenumber (nm or cm-1)
    REAL*8,INTENT(OUT)       :: wave_max      !Maximum wavelength or wavenumber (nm or cm-1)
    REAL*8,INTENT(OUT)       :: dwave         !Wavelength or wavenumber stepsize (nm or cm-1)
    
    CHARACTER*2,INTENT(OUT)  :: wave_type     !Indicates input type: wavelength (wl) or wavenumber (wn)
    
    ! Internal ------------------------------------------------------
    INTEGER, PARAMETER     :: iunin = 1
    INTEGER                :: i
    LOGICAL                :: chatty

    chatty = .false.

    !-----------------------------------------------------------------------
    ! Open the input file and read the first lines:
    !-----------------------------------------------------------------------
    OPEN(iunin,file=infile)
    DO i=1,4
       READ(iunin,*,err=899,end=898)
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the GRID parameters:
    !-----------------------------------------------------------------------
    READ(iunin,*,err=899,end=898) nx
    READ(iunin,*,err=899,end=898) ny
    READ(iunin,*,err=899,end=898) nz
    READ(iunin,*,err=899,end=898) x_dim
    READ(iunin,*,err=899,end=898) y_dim
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)

    nz = nz + 1 !add extra vacuum layer

    IF (chatty) PRINT *,'nx = ',nx
    IF (chatty) PRINT *,'ny = ',ny
    IF (chatty) PRINT *,'nz = ',nz
    IF (chatty) PRINT *,'x_dim = ',x_dim
    IF (chatty) PRINT *,'y_dim = ',y_dim

    DO i=1,4
       READ(iunin,*,err=899,end=898)
    ENDDO    

    !-----------------------------------------------------------------------
    ! Read the MONTE CARLO parameters:
    !-----------------------------------------------------------------------
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)

    DO i=1,4
       READ(iunin,*,err=899,end=898)
    ENDDO

    IF (chatty) PRINT *,'No MC parameters'

    !-----------------------------------------------------------------------
    ! Read the GENERAL parameters:
    !-----------------------------------------------------------------------
    READ(iunin,*,err=899,end=898) nmat
    READ(iunin,*,err=899,end=898) wave_type
    READ(iunin,*,err=899,end=898) wave_min
    READ(iunin,*,err=899,end=898) wave_max
    READ(iunin,*,err=899,end=898) dwave
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)
    READ(iunin,*,err=899,end=898)

    IF (chatty) PRINT *,'nmat = ',nmat
    IF (chatty) PRINT *,'wave_type = ',wave_type
    IF (chatty) PRINT *,'wave_min = ',wave_min
    IF (chatty) PRINT *,'wave_max = ',wave_max
    IF (chatty) PRINT *,'dwave = ',dwave    

    DO i=1,4
       READ(iunin,*,err=899,end=898)
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the AEROSOL parameters:
    !-----------------------------------------------------------------------
    READ(iunin,*,err=899,end=898) nscatt_types    

    IF (chatty) PRINT *,'nscatt_types = ',nscatt_types

    IF (nscatt_types.EQ.0) THEN
       READ(iunin,*,err=899,end=898)
       READ(iunin,*,err=899,end=898)
    ELSE
      DO i=1,nscatt_types
         READ(iunin,*,err=899,end=898)
         READ(iunin,*,err=899,end=898)
      ENDDO
    ENDIF

    READ(iunin,*,err=899,end=898) n_scat_angles

    IF (chatty) PRINT *,'nscatt_types = ',nscatt_types
    IF (chatty) PRINT *,'n_scat_angles = ',n_scat_angles

    DO i=1,4
       READ(iunin,*,err=899,end=898)
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the LAYER parameters:
    !-----------------------------------------------------------------------    
    READ(iunin,*,err=899,end=898) ngases
    READ(iunin,*,err=899,end=898) nclouds

    IF (chatty) PRINT *,'ngases = ',ngases
    IF (chatty) PRINT *,'nclouds = ',nclouds

    !-----------------------------------------------------------------------
    ! Check if the inputs are correct:
    !-----------------------------------------------------------------------
    IF (nmat.LT.1 .OR. nmat.GT.4 .OR. nmat.EQ.2) THEN
       WRITE(*,*) 'rdsize: nmat should be 1, 3 or 4'
       STOP
    ENDIF

    IF (ngases .NE. 1) THEN
       PRINT *,'ngases=',ngases
       PRINT *,'ERROR: ngases should be 1'
       PRINT *,'If you want to exclude gas absorption, put the vmrs equal to 0.0.'
       STOP
    ENDIF

    IF (nscatt_types .LT. 1) THEN
       PRINT *,'nscatt_types=',nscatt_types
       PRINT *,'ERROR: nscatt_types should be 1 or larger than 1'
       PRINT *,'If you want to exclude the droplets, put baer equal to 0.0.'
       STOP
    ENDIF    

    !-----------------------------------------------------------------------
    ! Close the input file:
    !-----------------------------------------------------------------------
    CLOSE(iunin)

    RETURN
898 WRITE(*,*) 'rdsize: unexpected end of file encountered'
    WRITE(*,*) 'Please make sure a (correct) copy of your input file is located in $MONKI_HOME/INPUT/'
    STOP
899 WRITE(*,*) 'rdsize: ERROR reading file'
    STOP

  END SUBROUTINE rdsize

!+-----------------------------------------------------------------
! SUBROUTINE rdinput
!+-----------------------------------------------------------------
! Read the remaining variables from the input file  
!+-----------------------------------------------------------------
  SUBROUTINE rdinput(infile,nx,ny,nz,nscatt_types,ngases,nclouds,       &
    &                vmrgas,pres,temp,alt,                              &
    &                read_bmsca_from_file,bmsca_file,depol_fixed,       &
    &                gas_col_type,fix_depol,read_bmabs_from_file,bmabs_file,&
    &                background,coefile,coelut,                         &
    &                cloud_input_type,cloud_abundance,                  &
    &                N_PHOTONS_PER_COLUMN,N_ROUNDS,                     &
    &                save_path,save_scats,                              &
    &                stop_converged,tol_error,detector_location,        &
    &                source_location,register_location)

    IMPLICIT NONE

    ! Input -------------------------------------------------------
    CHARACTER(len=500),INTENT(IN) :: infile       !Name of the input file
    INTEGER,INTENT(IN)            :: nx,ny,nz     !Number of grid cells in x, y and z directions (-)
    INTEGER,INTENT(IN)            :: nscatt_types !Number of cloud droplet types (-)
    INTEGER,INTENT(IN)            :: ngases       !Number of absorbing gases (-)
    INTEGER,INTENT(IN)            :: nclouds      !Maximum number of different clouds in a layer (-)  

    ! Output -------------------------------------------------------
    INTEGER,INTENT(OUT)       :: N_ROUNDS             !Number of simulation rounds used (-)
    INTEGER,INTENT(OUT)       :: gas_col_type         !Indicates if hydrostat. equil. is used (1) or exp. density decay (0)
    INTEGER,INTENT(OUT)       :: detector_location    !Location of the detector, TOA (1) or surface (0)
    INTEGER,INTENT(OUT)       :: register_location    !Location where the result is stored, TOA (1) or surface (0)
    INTEGER,INTENT(OUT)       :: source_location      !Location of the light source, TOA (1) or surface (0)

    INTEGER*8,INTENT(OUT)     :: N_PHOTONS_PER_COLUMN !Number of photons user per column and per simulation round (-)

    REAL*8,INTENT(OUT)        :: pres(nz+1)           !Pressure at the levels (Pa)
    REAL*8,INTENT(OUT)        :: temp(nz+1)           !Temperature at the levels (K)
    REAL*8,INTENT(OUT)        :: alt(nz+1)            !Altitude at the levels (m)
    REAL*8,INTENT(OUT)        :: vmrgas(ngases,nz+1)  !Volume mixing ratios at the levels (ppm)
    REAL*8,INTENT(OUT)        :: depol_fixed          !Fixed depolariation factor (-)
    REAL*8,INTENT(OUT)        :: tol_error            !Maximum allowed reflectane or transmittance standard deviation
    REAL*8,INTENT(OUT)        :: cloud_abundance(nscatt_types,nx,ny,nz) !COT (-) or number column density ()

    CHARACTER*30,INTENT(OUT)  :: coefile(nscatt_types) !Start name of file with aerosol single scattering expansion coefficients
    CHARACTER*30,INTENT(OUT)  :: coelut(nscatt_types) !Start name of file with LUT for aerosol single scattering matrix elements
    CHARACTER*500,INTENT(OUT) :: bmsca_file,bmabs_file !Filenames with bmsca en bmabs values
    CHARACTER*5,INTENT(OUT)   :: background            !Planet type for depol. factor and refr. index: Earth or Venus
    CHARACTER*4,INTENT(OUT)   :: cloud_input_type      !Cloud input type: COT (baer) or number column density (ncol)

    LOGICAL,INTENT(OUT)       :: save_path             !Indicates if AMFs are saved
    LOGICAL,INTENT(OUT)       :: save_scats            !Indicates if scattering contributions are saved
    LOGICAL,INTENT(OUT)       :: stop_converged        !Indicates if simulation stops after convergence
    LOGICAL,INTENT(OUT)       :: read_bmsca_from_file  !Indicates if bmsca are read from file
    LOGICAL,INTENT(OUT)       :: fix_depol             !Indicates if fixed depolarization factor is used from input file 
    LOGICAL,INTENT(OUT)       :: read_bmabs_from_file  !Indicates if bmabs are read from file

    ! Internal -----------------------------------------------------
    INTEGER                  :: i,j,k,igas,icloud,iostat
    INTEGER                  :: iaer(nclouds)
    INTEGER                  :: ixmin(nclouds),ixmax(nclouds)
    INTEGER                  :: iymin(nclouds),iymax(nclouds)
    INTEGER, PARAMETER       :: iunin = 1
    REAL*8                   :: xpres,xalt,cloud_val(nclouds)

    REAL*8,PARAMETER :: eps = 1.0D-16

    !-----------------------------------------------------------------------
    ! Open the input file and read the first lines:
    !-----------------------------------------------------------------------
    OPEN(iunin,file=infile)
    DO i=1,4
       READ(iunin,*,err=999,end=998)
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the GRID parameters:
    !-----------------------------------------------------------------------
    READ(iunin,*,err=999,end=998) 
    READ(iunin,*,err=999,end=998) 
    READ(iunin,*,err=999,end=998)
    READ(iunin,*,err=999,end=998) 
    READ(iunin,*,err=999,end=998)    
    READ(iunin,*,err=999,end=998) source_location
    READ(iunin,*,err=999,end=998) detector_location
    READ(iunin,*,err=999,end=998) register_location

    !PRINT*,'! detector_location=',detector_location
    !STOP

    DO i=1,4
       READ(iunin,*,err=999,end=998)
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the MONTE CARLO parameters:
    !-----------------------------------------------------------------------
    READ(iunin,*,err=999,end=998) N_PHOTONS_PER_COLUMN
    READ(iunin,*,err=999,end=998) N_ROUNDS
    READ(iunin,*,err=999,end=998) stop_converged
    READ(iunin,*,err=999,end=998) TOL_ERROR

    IF ((N_ROUNDS .LT. 30) .AND. stop_converged) THEN
       PRINT *,'Error in rdsize(): stop_converged = .true., but N_ROUNDS &
               & should be at least 30 for a good estimation of the reflectane &
               & standard deviation. Please increase N_ROUNDS or set stop_converged &
               & equal to .false. in the input file.'
       STOP
    ENDIF

    IF (N_ROUNDS .LT. 30) THEN
       PRINT *,'Warning in rdsize(): N_ROUNDS is smaller than 30, but should &
              & be at least 30 for a good estimation of the reflectane standard &
              & deviation. The standard deviation will now be set equal to 999. &
              & Proceed by pushing ENTER or increase N_ROUNDS.'
       READ(*,*)
    ENDIF

    !PRINT *,'N_PHOTONS_PER_COLUMN = ',N_PHOTONS_PER_COLUMN
    !PRINT *,'N_ROUNDS = ',N_ROUNDS
    !PRINT *,'TOL_ERROR = ',TOL_ERROR
    !PRINT *,'RADIANCE_MAX = ',RADIANCE_MAX

    DO i=1,4
       READ(iunin,*,err=999,end=998)
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the GENERAL parameters:
    !-----------------------------------------------------------------------
    READ (iunin,*,err=999,end=998)
    READ (iunin,*,err=999,end=998)
    READ (iunin,*,err=999,end=998) 
    READ (iunin,*,err=999,end=998) 
    READ (iunin,*,err=999,end=998) 
    READ (iunin,*,err=999,end=998) read_bmsca_from_file
    READ (iunin,'(A23)',err=999,end=998) bmsca_file    
    READ (iunin,*,err=999,end=998) gas_col_type    
    READ (iunin,*,err=999,end=998) background
    READ (iunin,*,err=999,end=998) fix_depol
    READ (iunin,*,err=999,end=998) depol_fixed    
    READ (iunin,*,err=999,end=998) read_bmabs_from_file
    READ (iunin,'(A23)',err=999,end=998) bmabs_file 
    READ (iunin,*,err=999,end=998) cloud_input_type  
    READ (iunin,*,err=999,end=998) save_path
    READ (iunin,*,err=999,end=998) save_scats

    DO i=1,4
       READ(iunin,*,err=999,end=998)
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the AEROSOL parameters:
    !-----------------------------------------------------------------------
    READ(iunin,*,err=999,end=998)

    IF (nscatt_types.EQ.0) THEN
       READ(iunin,*,err=999,end=998)
       READ(iunin,*,err=999,end=998)
    ELSE
      DO i=1,nscatt_types
         !READ(iunin,'(A18)') coefile(i)
         !READ(iunin,'(A18)') coelut(i)
         READ(iunin,'(A23)') coefile(i)
         READ(iunin,'(A23)') coelut(i)

      ENDDO
    ENDIF

    READ(iunin,*,err=999,end=998)

    DO i=1,4
       READ(iunin,*,err=999,end=998)
    ENDDO

    !-----------------------------------------------------------------------
    ! Read the LAYER parameters:
    !-----------------------------------------------------------------------
    READ(iunin,*,err=999,end=998)
    READ(iunin,*,err=999,end=998)

    DO i=1,3
       READ(iunin,*,err=999,end=998)
    ENDDO

    !Loop over levels (nlevels = nlayers+1 = nz+1)
    DO k=1,nz+1

       !PRINT*,'k=',k

       IF ( k .LE. nz) THEN

          READ(iunin,*,err=999,end=998) j,xalt,xpres,temp(k),(vmrgas(igas,k), igas=1,ngases), &
                 &  (cloud_val(icloud),iaer(icloud),                                       &
                 &  ixmin(icloud),ixmax(icloud),iymin(icloud),iymax(icloud),icloud=1,nclouds)

          !PRINT *,'j=',j,' xalt=',xalt,' km, xpres=',xpres, ' bar','temp(k)=',temp(k), ' K'
          !PRINT *,'vmrgas(:,k)=',vmrgas(:,k),', baer_val=',baer_val,', iaer(k)=',iaer(k)
          !PRINT *,'ixmin=',ixmin,', ixmax=',ixmax,', iymin=',iymin,', iymax=',iymax

          cloud_abundance(:,:,:,k) = 0.D0 !baer or ncol, depending on cloud_input_type

          DO icloud=1,nclouds
             
             IF (iaer(icloud).LT.0 .OR. iaer(icloud).GT.nscatt_types) iaer(icloud)=0

             IF (cloud_val(icloud).GT.eps) THEN

                cloud_abundance(iaer(icloud),ixmin(icloud):ixmax(icloud),iymin(icloud):iymax(icloud),k) = &
            &   cloud_abundance(iaer(icloud),ixmin(icloud):ixmax(icloud),iymin(icloud):iymax(icloud),k) + cloud_val(icloud)

                IF ((ixmin(icloud) .LT. 1) .OR. (ixmin(icloud) .GT. nx)) THEN
                    PRINT *,'Error: ixmin should not be smaller than 1 and not be larger than nx.'
                    STOP
                ELSEIF ((ixmax(icloud) .LT. 1) .OR. (ixmax(icloud) .GT. nx)) THEN
                   PRINT *,'Error: ixmax should not be smaller than 1 and not be larger than nx.'
                   STOP
                ELSEIF ((iymin(icloud) .LT. 1) .OR. (iymin(icloud) .GT. ny)) THEN
                   PRINT *,'Error: iymin should not be smaller than 1 and not be larger than ny.'
                STOP
                ELSEIF ((iymax(icloud) .LT. 1) .OR. (iymax(icloud) .GT. ny)) THEN
                   PRINT *,'Error: iymax should not be smaller than 1 and not be larger than ny.'
                   STOP
                ENDIF

             ENDIF

          ENDDO

          pres(k) = xpres*1.D5 !in Pa
          alt(k)  = xalt*1.D3  !in m        

          !Error
          IF (pres(k) .LT. eps) THEN
             PRINT *,'Error in rdinput: pressure equals zero.'
             PRINT *,'k=',k
             PRINT *,'Please provide a pressure larger than zero in the input file.'
             STOP
          ENDIF

          !Error
          IF (temp(k) .LT. eps) THEN
             PRINT *,'Error in rdinput: temperature equals zero.'
             PRINT *,'k=',k
             PRINT *,'Please provide a temperature larger than zero in the input file.'
             STOP
          ENDIF          

       ELSE

          !Add vacuum layer
          pres(k) = 0.D0
          alt(k)  = alt(k-1)  + 0.000001D0 !+ 300.D0*1.0D3
          temp(k) = temp(k-1) - 100.D0

       ENDIF

    ENDDO

    !PRINT *,'pres=',pres

    READ(iunin,*,IOSTAT=iostat) j
    IF (iostat .EQ. 0) THEN
       PRINT *,'Error in rdinput: there are more lines in the input file than levels (nz+1)'
       PRINT *,'Please make sure that nz = nlevels - 1 in the input file, and that' 
       PRINT *,'the input file ends after the last line.' 
       STOP
    ENDIF

    !-----------------------------------------------------------------------
    ! Close the input file:
    !-----------------------------------------------------------------------
    CLOSE(iunin)

    RETURN
998 WRITE(*,*) 'rdinput: unexpected end of file encountered'
    STOP
999 WRITE(*,*) 'rdinput: ERROR reading file'
    STOP
   
  END SUBROUTINE rdinput

END MODULE main_routines_cartesian

