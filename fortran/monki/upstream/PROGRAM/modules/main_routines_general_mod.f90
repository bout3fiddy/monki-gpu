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
!-------------------------------------------------------------------
! MODULE main_routines_general
!
! Module that contains subroutines for the main.f90 program of MONKI
! Subroutines that are independent of the grid type. 
!
! Author: Victor Trees
!-------------------------------------------------------------------
MODULE main_routines_general

  IMPLICIT NONE

CONTAINS

!+-----------------------------------------------------------------
! SUBROUTINE molecules_earth
!+-----------------------------------------------------------------
! Compute the optical properties of the gas
! This version is for the Earth's atmosphere
!+-----------------------------------------------------------------
  SUBROUTINE molecules_earth(nz,pres,temp,alt,wav,depol,gas_col_type,bmsca,monki_home,ngases,vmrgas,bmabs)

    USE physical_parameters
 
    IMPLICIT NONE

    ! Input -------------------------------------------------------
    INTEGER,INTENT(IN) :: nz            !Number of layers in the z direction (-)
    INTEGER,INTENT(IN) :: ngases        !Number of absorbing gases (-)
    INTEGER,INTENT(IN) :: gas_col_type  !Indicates if hydrostat. equil. if used (1) or exp. density decay (0)
    
    REAL*8,INTENT(IN)  :: pres(nz+1)    !Pressure at the levels (Pa)
    REAL*8,INTENT(IN)  :: wav           !Wavelength of the light (nm)
    REAL*8,INTENT(IN)  :: temp(nz+1)    !Temperature at the levels (K)
    REAL*8,INTENT(IN)  :: alt(nz+1)     !Altitude at the levels (m)
    REAL*8,INTENT(IN)  :: vmrgas(ngases,nz+1) !Volume mixing ratios of the absorbing gases at the levels (ppm)
    
    CHARACTER(len=300),INTENT(IN)  :: monki_home !MONKI home directory

    ! Output ------------------------------------------------------
    REAL*8,INTENT(OUT) :: depol     !Depolarization factor
    REAL*8,INTENT(OUT) :: bmsca(nz) !Gas scattering optical thicknesses in the layers (-)
    REAL*8,INTENT(OUT) :: bmabs(nz) !Gas absorption optical thicknesses in the layers (-)

    ! Internal ----------------------------------------------------
    INTEGER  :: k,igas
    REAL*8   :: nmol(nz),w2,ri,dep,rindex,ssca,wavmic
    REAL*8   :: pt(nz),pt1,pt2,dz
    REAL*8   :: bgasabs(ngases,nz),sigmu1,sigmu2
    REAL*8   :: sigmaabs(ngases,nz+1),ptsigmu(nz)

    REAL*8 :: King_factor_air,sig,wa
    REAL*8 :: KingFactor_N2,KingFactor_O2,KingFactor_Ar,KingFactor_CO2
    
    ! fractions are expressed as parts per volume by percent
    REAL*8,PARAMETER :: fraction_N2  =  78.084d0
    REAL*8,PARAMETER :: fraction_O2  =  20.9476d0
    REAL*8,PARAMETER :: fraction_Ar  =  0.934D0
    REAL*8,PARAMETER :: fraction_CO2 =  0.0375D0   ! for 375 ppmv CO2
    REAL*8,PARAMETER :: w_N2  =  28.0134D0 !g/mol
    REAL*8,PARAMETER :: w_O2  =  32.00D0   !g/mol
    REAL*8,PARAMETER :: w_Ar  =  39.948D0  !g/mol
    REAL*8,PARAMETER :: w_CO2 =  44.01D0   !g/mol
    REAL*8,PARAMETER :: grav  =  9.807D0   !m/s2 

    !REAL*8,PARAMETER :: eps = 1.D-10
    REAL*8,PARAMETER :: eps = 1.D-16

    CHARACTER*72 :: absfile

    LOGICAL  :: chatty,hydrostatic_equilibrium

    chatty = .false.

    !-----------------------------------------------------------------------
    ! Calculate the number of molecules in the layer (in m^-2):
    !-----------------------------------------------------------------------
    IF (gas_col_type .EQ. 1) THEN
       hydrostatic_equilibrium = .true.
    ELSEIF (gas_col_type .EQ. 0) THEN
       hydrostatic_equilibrium = .false.
    ELSE
       PRINT *,'ERROR: Wrong gas_col_type.'
       STOP
    ENDIF

    IF (hydrostatic_equilibrium) THEN

       !hydrostatic equilibrium
       !here we don't need the temperature profile
       !but we use the molecular mass

       !molecular mass in kg
       wa = (w_N2*fraction_N2 + w_O2*fraction_O2 + w_Ar*fraction_Ar + w_CO2*fraction_CO2)/ &
   &           (fraction_N2 + fraction_O2 + fraction_Ar + fraction_CO2)*1.D-3/avogad

       !loop over layers
       DO k=1,nz
          IF (k.EQ.nz) THEN
             nmol(k) = 0.D0 !vacuum layer
          ELSE
             !number density in 1/m2
             nmol(k)= (pres(k)-pres(k+1))/(wa*grav)
             IF (chatty) PRINT *, 'nmol(k)=',nmol(k),'pres(k)=',pres(k),'pres(k+1)=',pres(k+1)
          ENDIF
       ENDDO

    ELSE

       !no hydrostratic equilibrium (approach of DAK)
       !here we don't need the molecular mass calculation
       !but we use the given temperature profile

       !loop over layers
       DO k=1,nz

          IF (k.EQ.nz) THEN
             nmol(k) = 0.D0 !vacuum layer
          ELSE

             pt1 = pres(k)/temp(k)
             pt2 = pres(k+1)/temp(k+1)
             dz = alt(k+1)-alt(k)
             IF (DABS(pt1-pt2).LT.eps) THEN
                pt(k) = pt1*dz
             ELSE
                pt(k) = (pt2-pt1)*dz/(DLOG(pt2)-DLOG(pt1))
             ENDIF

             !number density in 1/m2
             nmol(k)= avogad/gasconst*pt(k)

          ENDIF
       ENDDO

    ENDIF

    !-----------------------------------------------------------------------
    ! Calculate the wavelength in micron instead of nanometer
    !-----------------------------------------------------------------------    
    wavmic = wav/1000.D0

    !-----------------------------------------------------------------------
    ! Calculate the refractive index "ri" for dry air at STP (288.15 K):
    ! Peck and Reader (1972)
    !
    ! NOTE that there is a discontinuity at wavel=0.16 microns:
    ! Just before this wavel, ri decreases steeply, and then increases
    ! sharply to gradually decrease again!
    !-----------------------------------------------------------------------
    IF (DABS(0.16D0-wavmic).LT.0.019D0) THEN
       WRITE(*,*)
       WRITE(*,*) 'WARNING: Discontinuity in air refractive index!'
       WRITE(*,*)
    ENDIF

    w2= 1.D0 /(wavmic*wavmic)
    !IF (wavmic.GT.0.23D0) THEN
    IF (.false.) THEN
       ri= 1.D0 +                                   &
      &    1.D-8*( 5791817.D0/ (238.0185D0 - w2)) + &
      &    1.D-8*( 167909.D0/ (57.362D0 - w2))
    ELSE
       ri= 1.D0 + 1.D-8*8060.51D0 +                 &
      &    1.D-8*( 2480990.D0/ (132.274D0 - w2) ) + &
      &    1.D-8*( 17455.7D0/ (39.32957D0 - w2) )
    ENDIF

    IF (chatty) PRINT *, 'ri=',ri    

    !ri = 1.00028698D0
    !PRINT *, 'use the one of Piet: ri=',ri

    !READ(*,*)

    !calculate the refractive index factor
    rindex= (ri*ri-1.D0)*(ri*ri-1.D0)/((ri*ri+2.D0)*(ri*ri+2.D0))
      
    !-----------------------------------------------------------------------
    ! Calculate the King factor
    ! 
    ! Based on the paper: "On Rayleigh Optical Depth Calculations" by
    ! BARRY A. BODHAINE, NORMAN B. WOOD, ELLSWORTH G. DUTTON and JAMES R. SLUSSER
    ! JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY, Vol.16, pp. 1854-1861, 1999.
    ! - no correction for water vapor is applied - therfore we have dry air here
    !-----------------------------------------------------------------------
    !dep= (6.D0+3.D0*depol)/(6.D0-7.D0*depol)

    sig = 1.0D3/wav

    KingFactor_N2  = 1.034d0 + 3.17d-4 * sig**2
    KingFactor_O2  = 1.096D0 + 1.385D-3 * sig**2 + 1.448D-4 * sig**4
    KingFactor_Ar  = 1.0D0
    KingFactor_CO2 = 1.15D0

    King_factor_air = fraction_N2  * KingFactor_N2  &
                    + fraction_O2  * KingFactor_O2  &
                    + fraction_Ar  * KingFactor_Ar  &
                    + fraction_CO2 * KingFactor_CO2

    dep = King_factor_air / ( fraction_N2 + fraction_O2 + fraction_Ar + fraction_CO2)

    IF (chatty) PRINT *, 'dep=',dep

    !Calculate depolarization factor
    depol = (6.D0*dep-6.D0)/(3.D0+7.D0*dep)

    IF (chatty) PRINT *, 'depol=',depol

    !-----------------------------------------------------------------------
    ! Calculate the molecular scattering cross section (in m^2):
    !-----------------------------------------------------------------------
    ssca= 24.D0*pi*pi*pi*rindex*dep/                      &
        &  (wavmic*wavmic*wavmic*wavmic*loschmidt*loschmidt*1.D-24)

    IF (chatty) PRINT *, 'ssca=',ssca

    !---------------------------------------------------------------------
    ! Calculate the molecular scattering optical thickness of
    ! each layer (incl vacuum layer):
    !---------------------------------------------------------------------
    DO k=1,nz
       bmsca(k)= ssca*nmol(k)
       IF (chatty) PRINT *, 'bmsca(k)=',bmsca(k)
    ENDDO

    !-----------------------------------------------------------------------
    ! Get the gas absorption cross-sections in the layers (depends on temp)
    ! in cm2/molec and convert to m2/molec
    !-----------------------------------------------------------------------
    IF (ngases .NE. 1) THEN
       PRINT *,'ngases=',ngases
       PRINT *,'ERROR: ngases should be 1'
       PRINT *,'If you want to exclude gas absorption, please set all vmrs equal to 0.'
       STOP
    ENDIF

    !--- ozone ---
    igas = 1 !index for O3, other gases not yet implemented
    absfile = '/DATA/abs_data/O3T_BP.dat'
    CALL get_O3_crabs(monki_home,absfile,wav,nz,temp,sigmaabs(igas,:))

    !PRINT *,'sigmaabs=',sigmaabs
    !READ(*,*)

    !--- conversion from cm2/molec to m2/molec
    sigmaabs = sigmaabs / 100.D0 / 100.D0

    !-----------------------------------------------------------------------
    ! Calculate the absorption optical thickness, bgasabs, of each gas
    ! in the layers
    !-----------------------------------------------------------------------
    DO igas=1,ngases

       DO k=1,nz

          IF (k.EQ.nz) THEN
             bgasabs(igas,k) = 0.D0 !vacuum layer
          ELSE

             pt1 = pres(k)/temp(k)
             pt2 = pres(k+1)/temp(k+1)
             dz = alt(k+1)-alt(k)
             sigmu1 = vmrgas(igas,k)/1000000.D0*sigmaabs(igas,k)
             sigmu2 = vmrgas(igas,k+1)/1000000.D0*sigmaabs(igas,k+1)

             IF (DABS(pt1-pt2).LT.eps) THEN
                !ptsigmu(k) = pt1*dz
                ptsigmu(k) = 0.5D0*pt1*(sigmu1+sigmu2)*dz
             ELSE
                !pt(k) = (pt2-pt1)*dz/(DLOG(pt2)-DLOG(pt1))
                ptsigmu(k) = (pt2-pt1)*dz*(sigmu1-sigmu2)/DLOG(pt2/pt1)**2 &
                      & + pt2*dz*sigmu2/DLOG(pt2/pt1)                   & 
                      & - pt1*dz*sigmu1/DLOG(pt2/pt1)                    
             ENDIF

             !PRINT *,'pt1=',pt1,', pt2=',pt2
             !PRINT *,'dz=',dz
             !PRINT *,'ptsigmu(k)=',ptsigmu(k)
             !PRINT *,'vmrgas(igas,k)=',vmrgas(igas,k)
             !PRINT *,'sigmaabs(igas,k)=',sigmaabs(igas,k)
             !READ(*,*)  

             bgasabs(igas,k)= avogad/gasconst*ptsigmu(k)

          ENDIF

       ENDDO

    ENDDO

    !-----------------------------------------------------------------------
    ! Adding of 'gasabs' for all gases to find 'bmabs' per layer
    !-----------------------------------------------------------------------
    DO k=1,nz
       bmabs(k)=0.D0
       IF (ngases.GE.1) THEN
          DO igas=1,ngases
             bmabs(k) = bmabs(k) + bgasabs(igas,k)
             IF (chatty) PRINT *, 'bmabs(k)=',bmabs(k)
          ENDDO
       ENDIF
    ENDDO


   RETURN


  END SUBROUTINE molecules_earth

!+-----------------------------------------------------------------
! SUBROUTINE molecules_venus
!+-----------------------------------------------------------------
! Compute the optical properties of the gas
! This version is for the Venus atmosphere
!+-----------------------------------------------------------------
  SUBROUTINE molecules_venus(nz,pres,temp,alt,wav,depol,gas_col_type,bmsca,monki_home,ngases,vmrgas,bmabs)

    USE physical_parameters

    IMPLICIT NONE

    ! Input -------------------------------------------------------
    INTEGER,INTENT(IN) :: nz            !Number of layers in the z direction (-)
    INTEGER,INTENT(IN) :: ngases        !Number of absorbing gases (-)
    INTEGER,INTENT(IN) :: gas_col_type  !Indicates if hydrostat. equil. if used (1) or exp. density decay (0)

    REAL*8,INTENT(IN)  :: pres(nz+1)    !Pressure at the levels (Pa)
    REAL*8,INTENT(IN)  :: wav           !Wavelength of the light (nm)
    REAL*8,INTENT(IN)  :: temp(nz+1)    !Temperature at the levels (K)
    REAL*8,INTENT(IN)  :: alt(nz+1)     !Altitude at the levels (m)
    REAL*8,INTENT(IN)  :: vmrgas(ngases,nz+1) !Volume mixing ratios of the absorbing gases at the levels (ppm)

    CHARACTER(len=300),INTENT(IN)  :: monki_home !MONKI home directory

    ! Output ------------------------------------------------------
    REAL*8,INTENT(OUT) :: depol     !Depolarization factor
    REAL*8,INTENT(OUT) :: bmsca(nz) !Gas scattering optical thicknesses in the layers (-)
    REAL*8,INTENT(OUT) :: bmabs(nz) !Gas absorption optical thicknesses in the layers (-)

    ! Internal ----------------------------------------------------
    INTEGER  :: k,igas
    
    REAL*8   :: nmol(nz),ri,dep,rindex,ssca,wavmic
    REAL*8   :: pt(nz),pt1,pt2,dz,wavnr2,wa
    REAL*8   :: t1,t2,t3,t4,t5
    REAL*8   :: bgasabs(ngases,nz),sigmu1,sigmu2
    REAL*8   :: sigmaabs(ngases,nz+1),ptsigmu(nz)

    ! fractions are expressed as parts per volume by percent
    REAL*8,PARAMETER :: fraction_N2  =   3.5D0
    REAL*8,PARAMETER :: fraction_CO2 =  96.5D0  
    REAL*8,PARAMETER :: w_N2  =  28.0134D0 !g/mol
    REAL*8,PARAMETER :: w_CO2 =  44.01D0 !g/mol
    REAL*8,PARAMETER :: grav = 8.87D0 !m/s2 
    REAL*8,PARAMETER :: avogad_mitchell=6.022169D23 !avogadro constant of Mitchell Yzer
    REAL*8,PARAMETER :: loschmidt_mitchell=2.5462880905D25 !value of Mitchell for 288.15 K and 1013 hPa (instead of 1013.25 hPa)
    !REAL*8,PARAMETER :: pi_mitchell=3.141592D0 !value of Mitchell

    !REAL*8,PARAMETER :: eps = 1.D-10
    REAL*8,PARAMETER :: eps = 1.D-16

    CHARACTER*72 :: absfile

    LOGICAL  :: chatty,hydrostatic_equilibrium

    chatty = .false.

    !-----------------------------------------------------------------------
    ! Calculate the number of molecules in the layer (in m^-2):
    !-----------------------------------------------------------------------
    IF (gas_col_type .EQ. 1) THEN
       hydrostatic_equilibrium = .true.
    ELSEIF (gas_col_type .EQ. 0) THEN
       hydrostatic_equilibrium = .false.
    ELSE 
       PRINT *,'ERROR: Wrong gas_col_type.'
       STOP
    ENDIF 

    IF (hydrostatic_equilibrium) THEN

       !hydrostatic equilibrium (approach of Mitchell Yzer)
       !here we don't need the temperature profile
       !but we use the molecular mass

       !molecular mass in kg
   !    wa = (w_CO2*fraction_CO2 + w_N2*fraction_N2)/ &
   !&           (fraction_CO2 + fraction_N2)*1.D-3/avogad
       wa = (w_CO2*fraction_CO2 + w_N2*fraction_N2)/ &
   &        (fraction_CO2 + fraction_N2)*1.D-3/avogad_mitchell           

       IF (chatty) PRINT *,'w_CO2=',w_CO2
       IF (chatty) PRINT *,'w_N2=',w_N2
       IF (chatty) PRINT *,'fraction_CO2=',fraction_CO2
       IF (chatty) PRINT *,'fraction_N2=',fraction_N2
       !IF (chatty) PRINT *,'(w_CO2*fraction_CO2 + w_N2*fraction_N2)/(fraction_CO2 + fraction_N2)=',(w_CO2*fraction_CO2 + w_N2*fraction_N2)/ &
   !&        (fraction_CO2 + fraction_N2)
       IF (chatty) PRINT *,'avogad_mitchell=',avogad_mitchell
       IF (chatty) PRINT *,'wa=',wa

       !loop over layers
       DO k=1,nz
          IF (k.EQ.nz) THEN
             nmol(k) = 0.D0 !vacuum layer
          ELSE
             !number density in 1/m2
             nmol(k)= (pres(k)-pres(k+1))/(wa*grav)
             IF (chatty) PRINT *, 'nmol(k)=',nmol(k),'pres(k)=',pres(k),'pres(k+1)=',pres(k+1)
          ENDIF
       ENDDO

    ELSE   

       !no hydrostratic equilibrium (approach of DAK)
       !here we don't need the molecular mass calculation
       !but we use the given temperature profile

       !loop over layers
       DO k=1,nz
          IF (k.EQ.nz) THEN
             nmol(k) = 0.D0 !vacuum layer
          ELSE
             pt1 = pres(k)/temp(k)
             pt2 = pres(k+1)/temp(k+1)
             dz = alt(k+1)-alt(k)
             IF (DABS(pt1-pt2).LT.eps) THEN
                pt(k) = pt1*dz
             ELSE
                pt(k) = (pt2-pt1)*dz/(DLOG(pt2)-DLOG(pt1))
             ENDIF

             !number density in 1/m2
             !nmol(k)= avogad/gasconst*pt(k)
             nmol(k)= avogad_mitchell/gasconst*pt(k)
          ENDIF
       ENDDO

    ENDIF  

    !-----------------------------------------------------------------------
    ! Calculate the wavelength in micron instead of nanometer
    !-----------------------------------------------------------------------    
    wavmic = wav/1000.D0

    !compare to Mitchell
    !wavmic = 2.4844597045480521D0 
    !wavmic = 1.3651839858457724D0

    !-----------------------------------------------------------------------
    ! Calculate the refractive index "ri" for CO2
    ! (Eq. 13 of Sneep and Ubachs, 2005)
    !
    ! Note that there is a type in Eq. 13 of Sneep and Ubachs (2005):
    ! 1e6 should be 1e3, that's why you find below the 1e-3 extra.
    !-----------------------------------------------------------------------

    !conversion from wavelength to wavenumber squared
    wavnr2 = 1.D0/(wavmic*wavmic*1.D-8)

    IF (chatty) PRINT *,'wavmic=',wavmic

    IF (chatty) PRINT *,'wavnr2=',wavnr2

    !calculate refractive index
    t1= 5799.25D0/((128908.9D0**2 - wavnr2))
    t2= 120.05D0/((89223.8D0**2 - wavnr2))
    t3= 5.3334D0/((75037.5D0**2 - wavnr2))
    t4= 4.3244D0/((67837.7D0**2 - wavnr2))
    t5= 0.1218145D-4/((2418.136D0**2 - wavnr2))
    ri= 1.D0 + 1.1427D6*(t1 + t2 + t3 + t4 + t5)*1.D-3    

    IF (chatty) PRINT *, 'ri=',ri

    !READ(*,*)

    !calculate the refractive index factor
    rindex= (ri*ri-1.D0)*(ri*ri-1.D0)/((ri*ri+2.D0)*(ri*ri+2.D0))

    IF (chatty) PRINT *,'rindex=',rindex

    !-----------------------------------------------------------------------
    ! Calculate the King factor for CO2
    ! (Eq. 14 of Sneep and Ubachs, 2005)
    !-----------------------------------------------------------------------
    !dep= (6.D0+3.D0*depol)/(6.D0-7.D0*depol)

    !Calculate King factor
    dep = 1.1364D0 + 25.3D-12*wavnr2

    IF (chatty) PRINT *, 'dep=',dep

    !Calculate depolarization factor
    depol = (6.D0*dep-6.D0)/(3.D0+7.D0*dep)

    IF (chatty) PRINT *, 'depol=',depol

    !-----------------------------------------------------------------------
    ! Calculate the molecular scattering cross section (in m^2):
    !-----------------------------------------------------------------------
    !ssca= 24.D0*pi*pi*pi*rindex*dep/                      &
   !&  (wavmic*wavmic*wavmic*wavmic*loschmidt*loschmidt*1.D-24)
    ssca= 24.D0*pi*pi*pi*rindex*dep/                      &
   &  (wavmic*wavmic*wavmic*wavmic*loschmidt_mitchell*loschmidt_mitchell*1.D-24)

    IF (chatty) PRINT *, 'pi=',pi
    IF (chatty) PRINT *, 'loschmidt_mitchell=',loschmidt_mitchell
    IF (chatty) PRINT *, 'ssca=',ssca

    !---------------------------------------------------------------------
    ! Calculate the molecular scattering optical thickness of
    ! each layer (incl vacuum layer):
    !---------------------------------------------------------------------
    DO k=1,nz
       bmsca(k)= ssca*nmol(k)
       IF (chatty) PRINT *, 'bmsca(k)=',bmsca(k)
    ENDDO

    !READ(*,*)

    !PRINT *,'bmsca=',bmsca
    !READ(*,*)

    !***********************************************************************
    ! Until here, change lines below later
    !***********************************************************************

    !-----------------------------------------------------------------------
    ! Get the gas absorption cross-sections in the layers (depends on temp)
    ! in cm2/molec and convert to m2/molec
    !-----------------------------------------------------------------------
    IF (ngases .NE. 1) THEN
       PRINT *,'ngases=',ngases
       PRINT *,'ERROR: ngases should be 1'
       PRINT *,'If you want to exclude gas absorption, please set all vmrs equal to 0.'
       STOP
    ENDIF

    !--- ozone ---
    igas = 1 !index for O3, other gases not yet implemented
    absfile = '/DATA/abs_data/O3T_BP.dat'
    CALL get_O3_crabs(monki_home,absfile,wav,nz,temp,sigmaabs(igas,:))

    !PRINT *,'sigmaabs=',sigmaabs
    !READ(*,*)

    !--- conversion from cm2/molec to m2/molec
    sigmaabs = sigmaabs / 100.D0 / 100.D0

    !-----------------------------------------------------------------------
    ! Calculate the absorption optical thickness, bgasabs, of each gas
    ! in the layers
    !-----------------------------------------------------------------------
    DO igas=1,ngases

       DO k=1,nz

          IF (k.EQ.nz) THEN
             bgasabs(igas,k) = 0.D0 !vacuum layer
          ELSE

             pt1 = pres(k)/temp(k)
             pt2 = pres(k+1)/temp(k+1)
             dz = alt(k+1)-alt(k)
             sigmu1 = vmrgas(igas,k)/1000000.D0*sigmaabs(igas,k)
             sigmu2 = vmrgas(igas,k+1)/1000000.D0*sigmaabs(igas,k+1)

             IF (DABS(pt1-pt2).LT.eps) THEN
                !ptsigmu(k) = pt1*dz
                ptsigmu(k) = 0.5D0*pt1*(sigmu1+sigmu2)*dz
             ELSE
                !pt(k) = (pt2-pt1)*dz/(DLOG(pt2)-DLOG(pt1))
                ptsigmu(k) = (pt2-pt1)*dz*(sigmu1-sigmu2)/DLOG(pt2/pt1)**2 &
                      & + pt2*dz*sigmu2/DLOG(pt2/pt1)                   &
                      & - pt1*dz*sigmu1/DLOG(pt2/pt1)
             ENDIF

             !PRINT *,'pt1=',pt1,', pt2=',pt2
             !PRINT *,'dz=',dz
             !PRINT *,'ptsigmu(k)=',ptsigmu(k)
             !PRINT *,'vmrgas(igas,k)=',vmrgas(igas,k)
             !PRINT *,'sigmaabs(igas,k)=',sigmaabs(igas,k)
             !READ(*,*)  

             bgasabs(igas,k)= avogad/gasconst*ptsigmu(k)

          ENDIF

       ENDDO

    ENDDO

    !-----------------------------------------------------------------------
    ! Adding of 'gasabs' for all gases to find 'bmabs' per layer
    !-----------------------------------------------------------------------
    DO k=1,nz
       bmabs(k)=0.D0
       IF (ngases.GE.1) THEN
          DO igas=1,ngases
             bmabs(k) = bmabs(k) + bgasabs(igas,k)
             IF (chatty) PRINT *, 'bmabs(k)=',bmabs(k)
          ENDDO
       ENDIF
    ENDDO


   RETURN


  END SUBROUTINE molecules_venus

!+-----------------------------------------------------------------
! SUBROUTINE aerosols
!+-----------------------------------------------------------------
! Compute the optical properties of the aerosols
!+-----------------------------------------------------------------
  SUBROUTINE aerosols(nmat,nscatt_types,new_scatt_info_pol,&
                    & monki_home,coelut,coefile,wav,cloud_input_type,&
                    & cloud_abundance,nx,ny,nz,baer)

    USE new_data_types
    USE physical_parameters

    IMPLICIT NONE

    ! Input -------------------------------------------------------
    INTEGER,INTENT(IN) :: nx,ny,nz     !Number of grid cells in the x, y and z directions (-)
    INTEGER,INTENT(IN) :: nmat         !Size of the Stokes vector: 1, 3 or 4 (-)
    INTEGER,INTENT(IN) :: nscatt_types !Number of cloud droplet types (-)

    REAL*8,INTENT(IN)  :: cloud_abundance(nscatt_types,nx,ny,nz) !COT (-) or number column density ()
    REAL*8,INTENT(IN)  :: wav                                    !Wavelength of the light (nm)
    
    CHARACTER(len=300),INTENT(IN)  :: monki_home               !MONKI home directory
    CHARACTER*4,INTENT(IN) :: cloud_input_type                 !Cloud input type: COT (baer) or number column density (ncol)
    CHARACTER*30,DIMENSION(nscatt_types),INTENT(IN) :: coelut  !Start name of file(s) with aerosol single scattering expansion coefficients
    CHARACTER*30,DIMENSION(nscatt_types),INTENT(IN) :: coefile !Start name of file(s) with LUT for aerosol single scattering matrix elements
    
    ! In- and output ----------------------------------------------
    TYPE(new_scatterer_info_pol),Dimension(nscatt_types),INTENT(INOUT) :: new_scatt_info_pol !Object with single scattering properties of the droplets
    
    ! Output ------------------------------------------------------
    REAL*8,INTENT(OUT) :: baer(nscatt_types,nx,ny,nz) !Aerosol total optical thicknesses in the cells (baersca + baerabs) (-)

    ! Internal ----------------------------------------------------
    INTEGER            :: is,ia,iunit,iunit_low,iunit_high,iline,ia_rev
    
    REAL*8             :: F11,F22,F33,F44,F12,F34,scat_angle,dum2,ssalb_aer
    REAL*8             :: F11_low,F22_low,F33_low,F44_low,F12_low,F34_low,ssalb_aer_low
    REAL*8             :: F11_high,F22_high,F33_high,F44_high,F12_high,F34_high,ssalb_aer_high
    REAL*8             :: CDF_bound,PDF_bound,PDF_bound_old,Q_ext,G,sigma
    REAL*8             :: Q_ext_low,G_low,Q_ext_high,G_high,dwav,wav_low,wav_high
    REAL*8             :: cos_theta_sca,cos_theta_sca_old
    
    CHARACTER(len=256) :: mie_filename,mie_file,mie_LUT_filename,mie_LUT_file
    CHARACTER(len=256) :: mie_filename_low,mie_file_low,mie_LUT_filename_low,mie_LUT_file_low
    CHARACTER(len=256) :: mie_filename_high,mie_file_high,mie_LUT_filename_high,mie_LUT_file_high
    CHARACTER(len=5)   :: dumm
    CHARACTER(len=6)   :: swav,swav_low,swav_high,dummy

    REAL*8,PARAMETER :: eps = 1.0D-16
    
    LOGICAL :: chatty,interp_mie_files

    chatty = .false.

    !---------------------------------------------------------------------
    ! Initilize cloud optical thickness
    !---------------------------------------------------------------------    
    baer(:,:,:,:) = 0.D0

    !---------------------------------------------------------------------
    ! Store wavelength in string
    !---------------------------------------------------------------------   
    IF (chatty) PRINT*,'wav=',wav

    interp_mie_files = .FALSE.

    IF (DMOD(wav,1.D0) .LT. eps) THEN
    
       IF ((10000D0.GT.wav).AND.(wav.GE.1000.D0)) THEN
          WRITE(dummy,'(F5.3)') wav/1000.D0
          swav= TRIM(ADJUSTL('0'//dummy))
       ELSEIF ((1000.D0.GT.wav).AND.(wav.GE.100.D0)) THEN
          WRITE(dumm,'(F4.3)') wav/1000.D0
          swav= TRIM(ADJUSTL('00'//dumm))
       ELSE
          PRINT*,'Error in aerosols: I have trouble coverting wav to a string.'
          STOP
       ENDIF

       IF (chatty) PRINT*,'swav=',swav

    ELSE

       !Interpolate between Mie files which should be provided at 
       !1 nm spectral resolution grid
       wav_low = REAL(INT(wav),8)
       wav_high = wav_low + 1.D0
       dwav = wav - wav_low
       interp_mie_files = .TRUE.

       IF (chatty) PRINT*,'wav_low=',wav_low
       IF (chatty) PRINT*,'wav_high=',wav_high
       IF (chatty) PRINT*,'dwav=',dwav

       IF ((10000D0.GT.wav_low).AND.(wav_low.GE.1000.D0)) THEN
          WRITE(dummy,'(F5.3)') wav_low/1000.D0
          swav_low= TRIM(ADJUSTL('0'//dummy))
       ELSEIF ((1000.D0.GT.wav_low).AND.(wav_low.GE.100.D0)) THEN
          WRITE(dumm,'(F4.3)') wav_low/1000.D0
          swav_low= TRIM(ADJUSTL('00'//dumm))
       ELSE
          PRINT*,'Error in aerosols: I have trouble coverting wav to a string.'
          STOP
       ENDIF

       IF ((10000D0.GT.wav_high).AND.(wav_high.GE.1000.D0)) THEN
          WRITE(dummy,'(F5.3)') wav_high/1000.D0
          swav_high= TRIM(ADJUSTL('0'//dummy))
       ELSEIF ((1000.D0.GT.wav_high).AND.(wav_high.GE.100.D0)) THEN
          WRITE(dumm,'(F4.3)') wav_high/1000.D0
          swav_high= TRIM(ADJUSTL('00'//dumm))
       ELSE
          PRINT*,'Error in aerosols: I have trouble coverting wav to a string.'
          STOP          
       ENDIF  

       IF (chatty) PRINT*,'dumm=',dumm
       IF (chatty) PRINT*,'swav_low=',swav_low
       IF (chatty) PRINT*,'swav_high=',swav_high

    ENDIF

    !---------------------------------------------------------------------
    ! Read the Mie scattering files
    !---------------------------------------------------------------------
    IF (chatty) PRINT *,'Read Mie scattering info'

    DO is=1,nscatt_types

          !find mie LUT file
          IF (.NOT. interp_mie_files) THEN
             mie_LUT_filename = '/DATA/Mie_data/'//TRIM(ADJUSTL(coelut(is)))//swav
             IF (chatty) PRINT*,'mie_LUT_filename=',mie_LUT_filename
             mie_LUT_file = TRIM(ADJUSTL(monki_home))//TRIM(ADJUSTL(mie_LUT_filename))
             iunit = 101
          ELSE
             mie_LUT_filename_low = '/DATA/Mie_data/'//TRIM(ADJUSTL(coelut(is)))//swav_low
             mie_LUT_filename_high = '/DATA/Mie_data/'//TRIM(ADJUSTL(coelut(is)))//swav_high
             IF (chatty) PRINT*,'mie_LUT_filename_low=',mie_LUT_filename_low
             IF (chatty) PRINT*,'mie_LUT_filename_high=',mie_LUT_filename_high
             mie_LUT_file_low = TRIM(ADJUSTL(monki_home))//TRIM(ADJUSTL(mie_LUT_filename_low))
             mie_LUT_file_high = TRIM(ADJUSTL(monki_home))//TRIM(ADJUSTL(mie_LUT_filename_high))
             iunit_low  = 102
             iunit_high = 103
          ENDIF

          !PRINT *,mie_LUT_file

          !open mie LUT file
          IF (.NOT. interp_mie_files) THEN
             OPEN(iunit,file=mie_LUT_file,status='old')
          ELSE
             OPEN(iunit_low,file=mie_LUT_file_low,status='old')
             OPEN(iunit_high,file=mie_LUT_file_high,status='old')
          ENDIF

          IF (chatty) PRINT *,'new_scatt_info_pol(is)%n_scat_angles=',new_scatt_info_pol(is)%n_scat_angles
          
          !loop over scattering angles
          DO ia=1,new_scatt_info_pol(is)%n_scat_angles

             !read scattering matrix elements
             !PRINT*,'read scattering matrix elements'

             IF (.NOT. interp_mie_files) THEN
                READ(iunit,*) scat_angle,F11,F22,F33,F44,F12,F34,dum2
             ELSE
                READ(iunit_low,*) scat_angle,F11_low,F22_low,F33_low,F44_low,F12_low,F34_low,dum2
                READ(iunit_high,*) scat_angle,F11_high,F22_high,F33_high,F44_high,F12_high,F34_high,dum2
                F11 = (1.D0-dwav)*F11_low + dwav*F11_high
                F22 = (1.D0-dwav)*F22_low + dwav*F22_high
                F33 = (1.D0-dwav)*F33_low + dwav*F33_high
                F44 = (1.D0-dwav)*F44_low + dwav*F44_high
                F12 = (1.D0-dwav)*F12_low + dwav*F12_high
                F34 = (1.D0-dwav)*F34_low + dwav*F34_high
             ENDIF

             !PRINT *,'scat_angle=',scat_angle
             !PRINT*,dum1,F11,F22,F33,F44,F12,F34,dum2
             !PRINT*,'iw=',iw,'ir=',ir,'ia=',ia

             !store scattering angle
             !let cos(scat_angle) range from -1 to 1 instead of 1 to -1
             ia_rev = new_scatt_info_pol(is)%n_scat_angles - ia + 1
             new_scatt_info_pol(is)%scat_angles(ia_rev)=DCOS(scat_angle*DegtoRad)

             !store scattering matrix elements
             new_scatt_info_pol(is)%F11(ia_rev)=F11
             new_scatt_info_pol(is)%F22(ia_rev)=F22
             new_scatt_info_pol(is)%F33(ia_rev)=F33
             new_scatt_info_pol(is)%F12(ia_rev)=F12
             new_scatt_info_pol(is)%F34(ia_rev)=F34
             new_scatt_info_pol(is)%F44(ia_rev)=F44

          ENDDO


          !Build CDF

          !loop over scattering angles
          DO ia=1,new_scatt_info_pol(is)%n_scat_angles          

             !scattering angle
             cos_theta_sca=new_scatt_info_pol(is)%scat_angles(ia)

             !get scattering matrix elements and normalize them
             F11=new_scatt_info_pol(is)%F11(ia)
             F12=new_scatt_info_pol(is)%F12(ia)

             !bound function
             IF (nmat.GT.1) THEN
                !PDF_bound = F11+DABS(F12)*DSQRT(2.D0)
                PDF_bound = F11+DABS(F12)
             ELSE
                PDF_bound = F11
             ENDIF

             !cumulative distribution function of bound function
             IF (ia.EQ.1) THEN
                CDF_bound = 0.D0
             ELSE                
                CDF_bound = CDF_bound + (cos_theta_sca-cos_theta_sca_old)*&
                                      & (PDF_bound+PDF_bound_old)/2.D0
             ENDIF 
   
             !pass numbers to the next angle
             cos_theta_sca_old = cos_theta_sca
             PDF_bound_old = PDF_bound

             !store CDF for high scattering angle resolution
             new_scatt_info_pol(is)%CDF_bound(ia)=CDF_bound

             !PRINT *,'ia = ',ia,', scat angle = ',DACOS(cos_theta_sca)/DegtoRad,', CDF_bound=',CDF_bound
             !PRINT *,'CDF_bound=',CDF_bound
             !PRINT *,CDF_bound

          ENDDO

         
          !Divide by total integral of cos(scat_angle) over -1 to 1 to get final CDF
          new_scatt_info_pol(is)%CDF_bound(:) = &
                  & new_scatt_info_pol(is)%CDF_bound(:) / CDF_bound

          !Close LUT files
          IF (.NOT. interp_mie_files) THEN
             CLOSE(iunit)
          ELSE
             CLOSE(iunit_low)
             CLOSE(iunit_high)
          ENDIF

          !Find mie file
          IF (.NOT. interp_mie_files) THEN
             mie_filename = '/DATA/Mie_data/'//TRIM(ADJUSTL(coefile(is)))//swav
             mie_file = TRIM(ADJUSTL(monki_home))//TRIM(ADJUSTL(mie_filename))
             iunit = 301
          ELSE
             mie_filename_low  = '/DATA/Mie_data/'//TRIM(ADJUSTL(coefile(is)))//swav_low
             mie_filename_high = '/DATA/Mie_data/'//TRIM(ADJUSTL(coefile(is)))//swav_high
             mie_file_low  = TRIM(ADJUSTL(monki_home))//TRIM(ADJUSTL(mie_filename_low))
             mie_file_high = TRIM(ADJUSTL(monki_home))//TRIM(ADJUSTL(mie_filename_high))
             iunit_low = 302
             iunit_high = 303
          ENDIF

          !Open mie file and read single scattering albedo of droplets
          IF (.NOT. interp_mie_files) THEN
             OPEN(iunit,file=mie_file,status='old')
             DO iline=1,7
                READ(iunit,*)
             ENDDO
             READ(iunit,'(32X,E25.14)') ssalb_aer
             READ(iunit,*)
             READ(iunit,'(32X,E25.14)') Q_ext
             READ(iunit,*)
             READ(iunit,'(32X,E25.14)') G
             CLOSE(iunit)

          ELSE
             OPEN(iunit_low,file=mie_file_low,status='old')
             DO iline=1,7
                READ(iunit_low,*)
             ENDDO
             READ(iunit_low,'(32X,E25.14)') ssalb_aer_low
             READ(iunit_low,*)
             READ(iunit_low,'(32X,E25.14)') Q_ext_low
             READ(iunit_low,*)
             READ(iunit_low,'(32X,E25.14)') G_low
             CLOSE(iunit_low)             
             
             OPEN(iunit_high,file=mie_file_high,status='old')
             DO iline=1,7
                READ(iunit_high,*)
             ENDDO
             READ(iunit_high,'(32X,E25.14)') ssalb_aer_high
             READ(iunit_high,*)
             READ(iunit_high,'(32X,E25.14)') Q_ext_high
             READ(iunit_high,*)
             READ(iunit_high,'(32X,E25.14)') G_high
             CLOSE(iunit_high)          
     
             ssalb_aer   = (1.D0-dwav)*ssalb_aer_low + dwav*ssalb_aer_high
             Q_ext       = (1.D0-dwav)*Q_ext_low + dwav*Q_ext_high
             G           = (1.D0-dwav)*G_low + dwav*G_high

          ENDIF

          !Store single scattering albedo
          new_scatt_info_pol(is)%ssalb_aer=ssalb_aer

       !---------------------------------------------------------------------
       ! Compute aerosol optical thickness
       !---------------------------------------------------------------------
       IF (cloud_input_type .EQ. 'ncol') THEN !ncol = particles per µm² in the layer
          sigma = Q_ext * G !extinction cross section in µm²
          baer(is,:,:,:) = cloud_abundance(is,:,:,:) * sigma
       ELSEIF (cloud_input_type .EQ. 'baer') THEN
          baer(is,:,:,:) = cloud_abundance(is,:,:,:)
       ELSE 
          PRINT*,'Error in aerosols: cloud_input_type should be ncol or baer'
          STOP
       ENDIF

    ENDDO


  END SUBROUTINE aerosols

!+----------------------------------------------------------------
! SUBROUTINE ADVANCE
!+----------------------------------------------------------------
! Advance the photon to the next cell and calculate DIST, the 
! distance it travelled to get there.  If the photon crosses a
! horizontal boundary of the model, inject it into the opposite
! side.
!+----------------------------------------------------------------  
  SUBROUTINE ADVANCE(NX, NY, NZ, I, J, K, X, Y, Z, &
     &     I_SIGN, J_SIGN, K_SIGN, A_ABS, B_ABS, C_ABS, &
     &     A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT)

    IMPLICIT NONE

    ! Input -------------------------------------------------------
    INTEGER,INTENT(IN)     :: NX,NY,NZ,I_SIGN,J_SIGN,K_SIGN
    REAL*8,INTENT(IN)      :: A_ABS,B_ABS,C_ABS
    REAL*8,INTENT(IN)      :: A_INV,B_INV,C_INV,XUNIT,YUNIT,ZUNIT(NZ)

    ! In- and output  ---------------------------------------------
    INTEGER,INTENT(INOUT)  :: I,J,K
    REAL*8,INTENT(INOUT)   :: X,Y,Z
    
    ! Internal ----------------------------------------------------
    REAL*8                 :: DIST,DIST_X,DIST_Y,DIST_Z

    !--------------------------------------------------------------
    ! Get distances to x, y and z boundaries along photon path 
    ! and compute the minimum
    !--------------------------------------------------------------
    DIST_X = X * A_INV
    DIST_Y = Y * B_INV
    DIST_Z = Z * C_INV
    DIST   = MIN(DIST_X, DIST_Y, DIST_Z)

    !--------------------------------------------------------------
    ! Determine intersected boundary and update photon direction
    !--------------------------------------------------------------
    IF (DIST_X .EQ. DIST) THEN ! the photon reaches the x boundary
       I = MOD(I + I_SIGN + NX - 1, NX) + 1
       X = XUNIT
       Y = DMAX1(0.D0, Y - B_ABS * DIST)
       Z = DMAX1(0.D0, Z - C_ABS * DIST)
    ELSE IF (DIST_Y .EQ. DIST) THEN ! the photon reaches the y boundary
       J = MOD(J + J_SIGN + NY - 1, NY) + 1
       X = DMAX1(0.D0, X - A_ABS * DIST)
       Y = YUNIT
       Z = DMAX1(0.D0, Z - C_ABS * DIST)
    ELSE                      ! the photon reaches the z boundary
       K = K + K_SIGN
       X = DMAX1(0.D0, X - A_ABS * DIST)
       Y = DMAX1(0.D0, Y - B_ABS * DIST)
       IF (K .GE. 1 .AND. K .LE. NZ) Z = ZUNIT(K)
    ENDIF
    
    RETURN

  END SUBROUTINE ADVANCE

!+-----------------------------------------------------------------
! SUBROUTINE flag_blocked_sensor_path
!+-----------------------------------------------------------------
! Flag the grid cells whose direct paths to the detector 
! are blocked by grid cells with aerosols
!+-----------------------------------------------------------------
  SUBROUTINE flag_blocked_sensor_path(nx,ny,nz,n_angles,nscatt_types,xunit,yunit,zunit,&
                                    & vartheta,phi,baer,k_detector,blocked_sensor_path)

    USE physical_parameters

    IMPLICIT NONE

    ! Input -------------------------------------------------------
    INTEGER,INTENT(IN) :: nx,ny,nz      !Number of grid cells in the x, y and z directions (-)
    INTEGER,INTENT(IN) :: n_angles      !Number of viewing zenith angles (-)
    INTEGER,INTENT(IN) :: nscatt_types  !Number of cloud droplet types (-)
    INTEGER,INTENT(IN) :: k_detector    !Location in the layers where detector is (-)
    
    REAL*8,INTENT(IN)  :: vartheta(n_angles) !Viewing zenith angles (deg)
    REAL*8,INTENT(IN)  :: phi                !Viewing azimuth angle (deg)
    REAL*8,INTENT(IN)  :: baer(nscatt_types,nx,ny,nz) !Aerosol total optical thicknesses in the cells (-)
    
    ! In- and output ----------------------------------------------
    INTEGER,INTENT(INOUT) :: blocked_sensor_path(nx,ny,nz,n_angles) !Flags of blocked paths to detector (-)

    ! Internal ----------------------------------------------------
    INTEGER  :: I_SIGN,J_SIGN,K_SIGN,I_photon,J_photon,K_photon
    INTEGER  :: I,J,K,ix,iy,iz,iangles,nix,niy,niz
    
    REAL*8   :: A_ABS,B_ABS,C_ABS,A,B,C
    REAL*8   :: A_INV,B_INV,C_INV,XUNIT,YUNIT,ZUNIT(NZ)
    REAL*8   :: X,Y,Z,ixfrac,iyfrac,izfrac

    REAL*8,PARAMETER :: eps = 1.0D-16

  !---------------------------------------------------------------------
  ! Flag grid cells for which the path to the sensor is blocked by a 
  ! cloud or aerosols. 
  !---------------------------------------------------------------------

  !Initialize
  blocked_sensor_path(:,:,:,:) = 0

  nix = 10
  niy = 10
  niz = 10
  IF (nx .EQ. 1) nix = 1
  IF (ny .EQ. 1) niy = 1
  IF (nz .EQ. 1) niz = 1


  DO iangles=1,n_angles

     a = DSIN(DegtoRad*vartheta(iangles)) * DCOS(DegtoRad*phi)
     b = DSIN(DegtoRad*vartheta(iangles)) * DSIN(DegtoRad*phi)
     c = DCOS(DegtoRad*vartheta(iangles))

     i_sign = INT(SIGN(1.D0, a))     !sign of a
     j_sign = INT(SIGN(1.D0, b))     !sign of b
     k_sign = INT(SIGN(1.D0, c))     !sign of c
     a_abs  = DABS(a)                !|a|
     b_abs  = DABS(b)                !|b|
     c_abs  = DABS(c)                !|c|
     a_inv  = 1.D0 / MAX(a_abs, eps) !1/|a|
     b_inv  = 1.D0 / MAX(b_abs, eps) !1/|b|
     c_inv  = 1.D0 / MAX(c_abs, eps) !1/|c|

     !First find the cells with aerosol
     DO k=1,nz
        !PRINT *,'k=',k
        DO j=1,ny
           !PRINT *,'j=',j
           DO i=1,nx
              !PRINT *,'i=',i

              IF (SUM(baer(:,i,j,k)) .GT. eps) THEN

                 !PRINT *,'i=',i,'j=',j,'k=',k    

                 !This grid cell contains aerosol
                 !blocked_sensor_path(i,j,k,iangles) = 1

                 !For all corners
                 !DO ix=0,1
                 !   X = ix*(XUNIT-eps)+(1-ix)*eps
                 !   DO iy=0,1
                 !      Y = iy*(YUNIT-eps)+(1-iy)*eps
                 !      DO iz=0,1
                 !         Z = iz*(ZUNIT(k)-eps)+(1-iz)*eps

                 DO ix=0,nix
                    ixfrac = ix/10.D0
                    X = ixfrac*(XUNIT-eps)+(1.D0-ixfrac)*eps
                    DO iy=0,niy
                       iyfrac = iy/10.D0
                       Y = iyfrac*(YUNIT-eps)+(1.D0-iyfrac)*eps
                       DO iz=0,niz
                          izfrac = iz/10.D0
                          Z = izfrac*(ZUNIT(k)-eps)+(1.D0-izfrac)*eps

                          !PRINT *,'a=',a,'b=',b,'c=',c
                          !PRINT *,'ix=',ix,'iy=',iy,'iz=',iz
                          !PRINT *,'X=',X,'Y=',Y,'Z=',Z

                          !Save aerosol location before calling ADVANCE
                          I_photon = I
                          J_photon = J
                          K_photon = K
                         
                          !DO WHILE ((K_photon .GE. 1 .AND. K_photon .LE. NZ).AND.(k_photon .NE. k_detector))
                          DO WHILE (K_photon .GE. 1 .AND. K_photon .LT. NZ) !k = nz is vacuum layer

                             !This grid cell is covered by aerosol as seen from the sensor
                             !PRINT *,'i_photon=',i_photon,'j_photon=',j_photon,'k_photon=',k_photon,'iangles=',iangles

                             blocked_sensor_path(i_photon,j_photon,k_photon,iangles) = 1

                             CALL ADVANCE(NX, NY, NZ, I_photon, J_photon, K_photon, X, Y, Z, &
                         &              -I_SIGN, -J_SIGN, -K_SIGN, A_ABS, B_ABS, C_ABS, &
                         &              A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT)

                             !PRINT *,'X=',X,'Y=',Y,'Z=',Z
                             !PRINT *,'-i_sign=',-i_sign,'-j_sign=',-j_sign,'-k_sign=',-k_sign
                             !PRINT *,'a_abs=',a_abs,'b_abs=',b_abs,'c_abs=',c_abs 
                             !PRINT *,'i_photon=',i_photon,'j_photon=',j_photon,'k_photon=',k_photon
                             !PRINT *,'---------'

                          ENDDO

                          !Get aerosol location back after calling ADVANCE
                          !I = I_BAK
                          !J = J_BAK
                          !K = K_BAK

                       ENDDO
                    ENDDO
                 ENDDO

              ENDIF

           ENDDO
        ENDDO
     ENDDO
  ENDDO

  END SUBROUTINE flag_blocked_sensor_path

END MODULE main_routines_general

