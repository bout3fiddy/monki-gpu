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
! MODULE mc_routines_general
!
! Module that contains subroutines for the mc.f90 program of MONKI
!
! Author: Victor Trees
!------------------------------------------------------------------
MODULE mc_routines_general

  IMPLICIT NONE

CONTAINS

!+----------------------------------------------------------------
! SUBROUTINE GROUND_LAMB2
!+----------------------------------------------------------------
! Sample the Lambertian surface reflectance angle
!+----------------------------------------------------------------         
SUBROUTINE GROUND_LAMB2(A,B,C,phi)

  USE physical_parameters

  IMPLICIT NONE

  ! Output -------------------------------------------------------
  REAL*8,INTENT(OUT) :: A, B, C !New photon direction (-)
  REAL*8,INTENT(OUT) :: phi     !New photon azimuth angle (rad)
  
  ! Internal -----------------------------------------------------
  REAL   :: P, U  
  REAL*8 :: theta
  
  !---------------------------------------------------------------
  ! Get uniformly distributed random numbers between 0 and 1
  !---------------------------------------------------------------
  CALL RANDOM_NUMBER(P)
  CALL RANDOM_NUMBER(U)

  !---------------------------------------------------------------
  ! Directly sample reflection zenith angle and azimuth angle
  !---------------------------------------------------------------
  theta = ACOS(SQRT(P))
  phi = 2.0D0 * pi * U

  !---------------------------------------------------------------
  ! Compute new photon direction
  !---------------------------------------------------------------
  A = DSIN(theta) * DCOS(phi)
  B = DSIN(theta) * DSIN(phi)
  C = DCOS(theta)

  RETURN
END SUBROUTINE GROUND_LAMB2


!+----------------------------------------------------------------
! SUBROUTINE ray_scat_matrix
!+----------------------------------------------------------------
! Get elements of the single scattering matrix for anisotropic
! Rayleigh scattering. See Hansen and Travis (1974).
!+----------------------------------------------------------------
SUBROUTINE ray_scat_matrix(nmat,cos_theta_sca,depol,F_sca)

  IMPLICIT NONE

  ! Input --------------------------------------------------------  
  INTEGER,INTENT(IN)  :: nmat !Size of the Stokes vector: 1, 3 or 4 (-)
  REAL*8,INTENT(IN)   :: depol !Depolarization factor (-)
  REAL*8,INTENT(IN)   :: cos_theta_sca !Cosine of scattering angle (-)

  ! Output -------------------------------------------------------  
  REAL*8,INTENT(OUT),DIMENSION(nmat,nmat) :: F_sca !Single Rayleigh scattering matrix (-)

  ! Internal -----------------------------------------------------  
  REAL*8              :: sin_theta_sca,Delta,Delta_prime
  REAL*8,PARAMETER    :: eps=1.0D-16

  !---------------------------------------------------------------  
  ! Get sine of scattering angle
  !---------------------------------------------------------------  
  IF ((1.D0-cos_theta_sca**2) .LT. eps) THEN
    sin_theta_sca = 0.D0
  ELSE
    sin_theta_sca=SQRT(1.D0-cos_theta_sca**2)
  ENDIF

  !---------------------------------------------------------------
  ! Compute depolarization terms
  !---------------------------------------------------------------
  Delta = (1.D0 - depol) / (1.D0 + depol/2.D0)
  Delta_prime = (1.D0 - 2.D0*depol)/(1.D0 - depol)

  !---------------------------------------------------------------
  ! Initialize scattering matrix
  !---------------------------------------------------------------
  F_sca = 0.D0

  !---------------------------------------------------------------
  ! Compute elements
  !---------------------------------------------------------------
  F_sca(1,1) = 3.D0/4.D0*(1.D0+cos_theta_sca**2)
  IF (nmat.GT.1) THEN
     F_sca(1,2) = -3.D0/4.D0*sin_theta_sca**2
     F_sca(2,1) = F_sca(1,2)
     F_sca(2,2) = 3.D0/4.D0*(1.D0+cos_theta_sca**2)
     F_sca(3,3) = 3.D0/2.D0*cos_theta_sca
     IF (nmat.GT.3) THEN
        F_sca(4,4) = Delta_prime * F_sca(3,3)
     ENDIF
  ENDIF

  !---------------------------------------------------------------
  ! Multiplications by depolarization term
  !---------------------------------------------------------------
  F_sca = Delta * F_sca
  F_sca(1,1) = F_sca(1,1) + (1.D0 - Delta)

  !---------------------------------------------------------------
  ! Normalize scattering matrix
  ! (not divided by pi, such that here you actually multiply radiance by pi (flux))
  !---------------------------------------------------------------
  F_sca = F_sca / 4.D0

END SUBROUTINE ray_scat_matrix        

!+----------------------------------------------------------------
! SUBROUTINE ray_scat_elements
!+----------------------------------------------------------------
! Get elements of the single scattering matrix for anisotropic
! Rayleigh scattering. See Hansen and Travis (1974).
! Only elements F11 and F12.
!
! Victor Trees, December 2022
!+----------------------------------------------------------------
SUBROUTINE ray_scat_elements(cos_theta_sca,depol,F11,F12)

  IMPLICIT NONE

  ! Input --------------------------------------------------------
  REAL*8,INTENT(IN)   :: depol !Depolarization factor (-)
  REAL*8,INTENT(IN)   :: cos_theta_sca !Cosine of scattering angle (-)

  ! Output -------------------------------------------------------- 
  REAL*8,INTENT(OUT)  :: F11,F12 !Single Rayleigh scattering matrix elements (-)

  ! Internal --------------------------------------------------------
  REAL*8              :: sin_theta_sca,Delta
  REAL*8,PARAMETER    :: eps=1.0D-16

  !---------------------------------------------------------------
  ! Get sine of scattering angle
  !---------------------------------------------------------------
  IF ((1.D0-cos_theta_sca**2) .LT. eps) THEN
    sin_theta_sca = 0.D0
  ELSE
    sin_theta_sca=SQRT(1.D0-cos_theta_sca**2)
  ENDIF

  !---------------------------------------------------------------
  ! Compute depolarization term
  !---------------------------------------------------------------
  Delta = (1.D0 - depol) / (1.D0 + depol/2.D0)

  !---------------------------------------------------------------
  ! Compute elements
  !---------------------------------------------------------------
  F11 = 3.D0/4.D0*(1.D0+cos_theta_sca**2)
  F12 = -3.D0/4.D0*sin_theta_sca**2

  !---------------------------------------------------------------
  ! Multiplications by depolarization term
  !---------------------------------------------------------------
  F11 = Delta * F11
  F12 = Delta * F12
  
  F11 = F11 + (1.D0 - Delta)

  !---------------------------------------------------------------
  ! Normalize scattering matrix elements
  ! (not divided by pi, such that here you actually multiply radiance by pi (flux))
  !---------------------------------------------------------------
  F11 = F11 / 4.D0
  F12 = F12 / 4.D0

END SUBROUTINE ray_scat_elements

!+----------------------------------------------------------------
! SUBROUTINE ray_scat_element
!+----------------------------------------------------------------
! Get element of the single scattering matrix for anisotropic
! Rayleigh scattering. See Hansen and Travis (1974).
! Only element F11.
!
! Victor Trees, December 2024
!+----------------------------------------------------------------
SUBROUTINE ray_scat_element(cos_theta_sca,depol,F11)

  IMPLICIT NONE

  ! Input --------------------------------------------------------
  REAL*8,INTENT(IN)   :: depol !Depolarization factor (-)
  REAL*8,INTENT(IN)   :: cos_theta_sca !Cosine of scattering angle (-)

  ! Output --------------------------------------------------------
  REAL*8,INTENT(OUT)  :: F11 !Rayleigh scattering phase function (-)

  ! Internal ------------------------------------------------------
  REAL*8              :: Delta
  REAL*8,PARAMETER    :: eps=1.0D-16

  !---------------------------------------------------------------
  ! Compute depolarization term
  !---------------------------------------------------------------
  Delta = (1.D0 - depol) / (1.D0 + depol/2.D0)

  !---------------------------------------------------------------
  ! Compute elements
  !---------------------------------------------------------------
  F11 = 3.D0/4.D0*(1.D0+cos_theta_sca**2)

  !---------------------------------------------------------------
  ! Multiplications by depolarization term
  !---------------------------------------------------------------
  F11 = Delta * F11

  F11 = F11 + (1.D0 - Delta)

  !---------------------------------------------------------------
  ! Normalize scattering matrix elements
  ! (not divided by pi, such that here you actually multiply radiance by pi (flux))
  !---------------------------------------------------------------
  F11 = F11 / 4.D0

END SUBROUTINE ray_scat_element

!+-----------------------------------------------------------------
! SUBROUTINE mie_scat_matrix
!+----------------------------------------------------------------
! Get elements of the single scattering matrix for Mie scattering, 
! generated using the expansion coefficients of 
! de Rooij and van der Stap (1984).
!
! Victor Trees, September 2022
!+----------------------------------------------------------------
SUBROUTINE mie_scat_matrix(nmat,cos_theta_sca,new_scatt_data_pol,n_scat_angles,F_sca)

  USE physical_parameters
  USE new_data_types

  IMPLICIT NONE

  ! Input --------------------------------------------------------
  INTEGER,INTENT(IN) :: nmat !Size of the Stokes vector: 1, 3 or 4 (-)
  INTEGER,INTENT(IN) :: n_scat_angles !Number of single scattering angles (-)
  REAL*8,INTENT(IN)  :: cos_theta_sca !Cosine of the single scattering angle (-)
  TYPE(new_scatterer_info_pol),INTENT(IN) :: new_scatt_data_pol !Object with droplet single scattering info

  ! Output -------------------------------------------------------
  REAL*8,INTENT(OUT),DIMENSION(nmat,nmat) :: F_sca !Single Mie scattering matrix

  ! Internal -----------------------------------------------------
  INTEGER             :: imin,imax
  REAL*8              :: scat_angles(2)
  REAL*8              :: F11_angles(2),F22_angles(2),F33_angles(2)
  REAL*8              :: F44_angles(2),F12_angles(2),F34_angles(2)
  REAL*8              :: F11,F22,F33,F44,F12,F34
  REAL*8,PARAMETER    :: eps=1.0D-16

  !---------------------------------------------------------------
  ! Initialize scattering matrix
  !---------------------------------------------------------------
  F_sca = 0.D0

  !---------------------------------------------------------------
  ! Estimate the index of the minimum scattering angle
  !---------------------------------------------------------------
  imin = MAX(1,INT((180.D0 - DACOS(cos_theta_sca)/DegtoRad)/180.D0*(n_scat_angles-1))+1)
  imax = MIN(imin+1,n_scat_angles)

  !---------------------------------------------------------------
  IF (imax.EQ.imin) THEN
     F11 = new_scatt_data_pol%F11(imax)
     F_sca(1,1) = F11 !<<<<<<<<<<<<<<<<<<<
     IF (nmat.GT.1) THEN
        F22 = new_scatt_data_pol%F22(imax)
        F33 = new_scatt_data_pol%F33(imax)
        F12 = new_scatt_data_pol%F12(imax)
        F_sca(2,2) = F22 !<<<<<<<<<<<<<<<<<<<
        F_sca(3,3) = F33 !<<<<<<<<<<<<<<<<<<<
        F_sca(1,2) = F12 !<<<<<<<<<<<<<<<<<<<
        F_sca(2,1) = F12 !<<<<<<<<<<<<<<<<<<<
        IF (nmat.GT.3) THEN
           F34 = new_scatt_data_pol%F34(imax)
           F44 = new_scatt_data_pol%F44(imax)
           F_sca(3,4) = F34  !<<<<<<<<<<<<<<<<<<<
           F_sca(4,3) = -F34 !<<<<<<<<<<<<<<<<<<<
           F_sca(4,4) = F44  !<<<<<<<<<<<<<<<<<<<
        ENDIF
     ENDIF
  ELSE
          
     !Get single scattering angles
     scat_angles = new_scatt_data_pol%scat_angles(imin:imax)

     !Sometimes, cos_theta_sca can differ slightly (~16th decimal) from the
     !tabulated value as a result of the reverse calculation of the single scat angle above
     IF (cos_theta_sca .LT. scat_angles(1)) THEN
        !PRINT *,'Warning in mie_scat_matrix: cos_theta_sca .LT. scat_angles(1)'
        !PRINT *,'cos_theta_sca .LT. scat_angles(1)'
        !PRINT *,'cos_theta_sca =',cos_theta_sca
        !PRINT *,'scat_angles=',scat_angles
        imin = imin - 1
        imax = imax - 1
        scat_angles = new_scatt_data_pol%scat_angles(imin:imax)
     ELSE IF (cos_theta_sca .GT. scat_angles(2)) THEN
        !PRINT *,'Warning in mie_scat_matrix: cos_theta_sca .GT. scat_angles(2)'
        !PRINT *,'cos_theta_sca .LT. scat_angles(1)'
        !PRINT *,'cos_theta_sca =',cos_theta_sca
        !PRINT *,'scat_angles=',scat_angles
        imin = imin + 1
        imax = imax + 1
        scat_angles = new_scatt_data_pol%scat_angles(imin:imax)
     ENDIF

     IF (cos_theta_sca .LT. scat_angles(1)) THEN
        PRINT *,'Error in mie_scat_matrix!'
        PRINT *,'cos_theta_sca .LT. scat_angles(1)'
        PRINT *,'cos_theta_sca =',cos_theta_sca
        PRINT *,'scat_angles=',scat_angles
        STOP
     ENDIF
     IF (cos_theta_sca .GT. scat_angles(2)) THEN
        PRINT *,'Error in mie_scat_matrix!'
        PRINT *,'cos_theta_sca .GT. scat_angles(2)'
        PRINT *,'cos_theta_sca=',cos_theta_sca
        PRINT *,'scat_angles=',scat_angles
        STOP
     ENDIF

     !Get elements versus scattering angle
     
     F11_angles = new_scatt_data_pol%F11(imin:imax)
     CALL linpol_dp_single(scat_angles, F11_angles, 2, cos_theta_sca, F11)
     F_sca(1,1) = F11
     
     IF (nmat .GT.1) THEN
        F22_angles = new_scatt_data_pol%F22(imin:imax)
        F33_angles = new_scatt_data_pol%F33(imin:imax)
        F12_angles = new_scatt_data_pol%F12(imin:imax)
        CALL linpol_dp_single(scat_angles, F22_angles, 2, cos_theta_sca, F22)
        CALL linpol_dp_single(scat_angles, F33_angles, 2, cos_theta_sca, F33)
        CALL linpol_dp_single(scat_angles, F12_angles, 2, cos_theta_sca, F12)
        F_sca(2,2) = F22
        F_sca(3,3) = F33
        F_sca(1,2) = F12
        F_sca(2,1) = F12
     
        IF (nmat .GT. 3) THEN
           F34_angles = new_scatt_data_pol%F34(imin:imax)
           F44_angles = new_scatt_data_pol%F44(imin:imax)
           CALL linpol_dp_single(scat_angles, F34_angles, 2, cos_theta_sca, F34)
           CALL linpol_dp_single(scat_angles, F44_angles, 2, cos_theta_sca, F44)
           F_sca(3,4) = F34
           F_sca(4,3) = -F34
           F_sca(4,4) = F44
           
        ENDIF
     ENDIF

  ENDIF

  !Normalize scattering matrix
  !(not divided by pi, such that here you actually multiply radiance by pi (flux))
  F_sca = F_sca / 4.D0

END SUBROUTINE mie_scat_matrix

!+-----------------------------------------------------------------
! SUBROUTINE mie_scat_elements
!+-----------------------------------------------------------------
! Get elements of the single scattering matrix for Mie scattering,
! generated using the expansion coefficients of
! de Rooij and van der Stap (1984).
! Only elements F11 and F12.
!
! Victor Trees, December 2022
!+-----------------------------------------------------------------
SUBROUTINE mie_scat_elements(cos_theta_sca,new_scatt_data_pol,n_scat_angles,F11,F12)

  USE physical_parameters
  USE new_data_types

  IMPLICIT NONE

  ! Input --------------------------------------------------------
  INTEGER,INTENT(IN) :: n_scat_angles !Number of single scattering angles (-)
  REAL*8,INTENT(IN)   :: cos_theta_sca !Cosine of the single scattering angle (-)
  TYPE(new_scatterer_info_pol),INTENT(IN) :: new_scatt_data_pol !Object with droplet single scattering info

  ! Output --------------------------------------------------------
  REAL*8,INTENT(OUT)  :: F11,F12 !Single Mie scattering matrix elements (-)

  ! Internal ------------------------------------------------------
  INTEGER             :: imin,imax
  REAL*8              :: F11_angles(2),F12_angles(2)
  REAL*8              :: scat_angles(2)  
  REAL*8,PARAMETER    :: eps=1.0D-16


  !estimate the index of the minimum scattering angle
  imin = MAX(1,INT((180.D0 - DACOS(cos_theta_sca)/DegtoRad)/180.D0*(n_scat_angles-1))+1)
  imax = MIN(imin+1,n_scat_angles)

  IF (imax.EQ.imin) THEN
     !PRINT *,'imax=imin'
     !PRINT *,'imax=',imax,'imin=',imin
     !READ(*,*)
     F11 = new_scatt_data_pol%F11(imax)
     F12 = new_scatt_data_pol%F12(imax)
  ELSE

     !Get single scattering angles
     scat_angles = new_scatt_data_pol%scat_angles(imin:imax)

     !Sometimes, cos_theta_sca can differ slightly (~16th decimal) from the
     !tabulated value as a result of the reverse calculation of the single scat angle above     
     IF (cos_theta_sca .LT. scat_angles(1)) THEN
        !PRINT *,'Warning in mie_scat_elements: cos_theta_sca .LT. scat_angles(1)'
        imin = imin - 1
        imax = imax - 1
        scat_angles = new_scatt_data_pol%scat_angles(imin:imax)
     ELSE IF (cos_theta_sca .GT. scat_angles(2)) THEN
        !PRINT *,'Warning in mie_scat_elements: cos_theta_sca .GT. scat_angles(2)'
        imin = imin + 1
        imax = imax + 1
        scat_angles = new_scatt_data_pol%scat_angles(imin:imax)
     ENDIF

     IF (cos_theta_sca .LT. scat_angles(1)) THEN
        PRINT *,'Error in mie_scat_elements!'
        PRINT *,'cos_theta_sca .LT. scat_angles(1)'
        PRINT *,'cos_theta_sca =',cos_theta_sca
        PRINT *,'scat_angles=',scat_angles
        STOP
     ENDIF
     IF (cos_theta_sca .GT. scat_angles(2)) THEN
        PRINT *,'Error in mie_scat_elements!'
        PRINT *,'cos_theta_sca .GT. scat_angles(2)'
        PRINT *,'cos_theta_sca=',cos_theta_sca
        PRINT *,'scat_angles=',scat_angles
        STOP
     ENDIF

     !Get elements versus scattering angle
     F11_angles = new_scatt_data_pol%F11(imin:imax)
     F12_angles = new_scatt_data_pol%F12(imin:imax)

     CALL linpol_dp_single(scat_angles, F11_angles, 2, cos_theta_sca, F11)
     CALL linpol_dp_single(scat_angles, F12_angles, 2, cos_theta_sca, F12)

  ENDIF

  !Normalize scattering matrix elements
  !(not divided by pi, such that here you actually multiply radiance by pi (flux))
  F11 = F11 / 4.D0
  F12 = F12 / 4.D0

END SUBROUTINE mie_scat_elements

!+-----------------------------------------------------------------
! SUBROUTINE rayl_scatter_pol
!+-----------------------------------------------------------------
! Sample the cosine of a photon scattering angle and rotation angle
! using rejection sampling,
! and compute scattering zenith angle, azimuth difference,
! and second rotation angle, for anisotropic Rayleigh scattering.
!
! Victor Trees, September 2022
!+-----------------------------------------------------------------
SUBROUTINE rayl_scatter_pol(backward_mode,nmat,depol,S_in,cos_theta_in,cos_theta_out,dphi,cos_theta_sca)

  !USE knuth_random

  IMPLICIT NONE

  ! Input --------------------------------------------------------
  INTEGER,INTENT(IN) :: nmat !Size of the Stokes vector: 1, 3 or 4 (-)
  REAL*8,INTENT(IN)  :: depol !Depolarization factor (-)
  REAL*8,INTENT(IN)  :: cos_theta_in !Cosine of incident zenith angle (-)
  REAL*8,INTENT(IN)  :: S_in(nmat) !Incident Stokes vector (-)
  LOGICAL            :: backward_mode !Indicates if we run in backward_mode

  ! Output -------------------------------------------------------
  REAL*8,INTENT(OUT) :: cos_theta_out !Cosine of scattered zenith angle (-)
  REAL*8,INTENT(OUT) :: dphi !Azimuth difference between incident and scattered beams (rad)
  REAL*8,INTENT(OUT) :: cos_theta_sca !Cosine of single scattering angle (-)

  ! Internal -----------------------------------------------------
  REAL*8 :: sigma1 !polarization rotation angle
  REAL*8 :: PDF,PDF_max,point,frac
  REAL   :: ran1,ran2,ran3
  REAL*8 :: ran_cos_theta_sca,ran_sigma1
  REAL*8 :: sin_theta_in,sin_theta_out,sin_theta_sca
  REAL*8 :: F11,F12,sign

  REAL*8,PARAMETER :: pi=3.141592653589793D0
  REAL*8,PARAMETER :: EPS=1.0D-16

  LOGICAL             :: chatty,bruteforce

  chatty = .false.
  bruteforce = .false.

  !-----------------------------------------------------------------------
  ! 1.) First calculate the maximum PDF
  !----------------------------------------------------------------------- 
  IF (chatty) PRINT *,'First calculate the maximum PDF'

  !-----------------------------------------------------------------------
  ! Get first Rayleigh scattering matrix element for scat angle of 0 deg
  ! Note: we assume that the PDF is max here (we have safety check later)
  !----------------------------------------------------------------------- 
  CALL ray_scat_element(1.D0,depol,F11)

  IF (chatty) PRINT *,'F11=',F11
  
  !-----------------------------------------------------------------------
  ! Brute force calculate sigma1 for which PDF is max and save PDF_max
  !-----------------------------------------------------------------------
  PDF_max = F11
  
  IF (chatty) PRINT *,'PDF_max=',PDF_max

  !-----------------------------------------------------------------------
  ! 2.) Secondly, sample cos_theta_sca and sigma1
  !----------------------------------------------------------------------- 
  IF (chatty) PRINT *,' Secondly, sample cos_theta_sca and sigma1'

  !-----------------------------------------------------------------------
  ! Sample random cosine of scattering angle between -1 and 1 
  !----------------------------------------------------------------------- 
100  CALL RANDOM_NUMBER(ran1)  
  !CALL ran(ran1)
  ran_cos_theta_sca = 2.D0*ran1 - 1.D0

  IF (chatty) PRINT *,'ran_cos_theta_sca',ran_cos_theta_sca

  !-----------------------------------------------------------------------
  ! Sample random rotation angle sigma1 between -pi and pi
  !-----------------------------------------------------------------------   
  !CALL ran(ran2)
  CALL RANDOM_NUMBER(ran2)
  ran_sigma1 = 2.D0*pi*ran2 - pi 

  IF (chatty) PRINT *,'ran_sigma1',ran_sigma1

  !-----------------------------------------------------------------------
  ! Get Rayleigh scattering matrix for this scattering angle
  !----------------------------------------------------------------------- 
  !CALL ray_scat_matrix(ran_cos_theta_sca,depol,F_sca) 
  CALL ray_scat_elements(ran_cos_theta_sca,depol,F11,F12)  
  IF (nmat .EQ. 1) F12 = 0.D0

  !IF (chatty) PRINT *,'F_sca',F_sca
  IF (chatty) PRINT *,'F11',F11
  IF (chatty) PRINT *,'F12',F12

  !-----------------------------------------------------------------------
  ! Set sign
  !-----------------------------------------------------------------------
  IF (backward_mode) THEN
     sign = 1.D0
  ELSE
     sign = -1.D0
  ENDIF

  !-----------------------------------------------------------------------
  ! Compute PDF(costheta,sigma1) for anisotropic Rayleigh scattering
  !-----------------------------------------------------------------------  
  IF (nmat .EQ. 1) THEN
     !PDF = F11*S_in(1)
     PDF = F11
  ELSE
     !PDF = F11*S_in(1) + F12*(S_in(2)*DCOS(2.D0*ran_sigma1)-S_in(3)*DSIN(2.D0*ran_sigma1))
     !PDF = F11 + F12*(S_in(2)*DCOS(2.D0*ran_sigma1)-S_in(3)*DSIN(2.D0*ran_sigma1))
     PDF = F11 + F12*(S_in(2)*DCOS(2.D0*ran_sigma1)+sign*S_in(3)*DSIN(2.D0*ran_sigma1))
  ENDIF     

  IF (chatty) PRINT *,'PDF',PDF

  !-----------------------------------------------------------------------
  ! Generate random point between 0 and max of PDF
  !-----------------------------------------------------------------------    
  !CALL ran(ran3)
  CALL RANDOM_NUMBER(ran3)
  point = ran3*PDF_max

  IF (chatty) PRINT *,'point',point

  !-----------------------------------------------------------------------
  ! Reject cos_theta_sca and sigma1 if point is larger than PDF
  !-----------------------------------------------------------------------   
  IF (point .GT. PDF) THEN

    IF (chatty) PRINT *,'Reject. GOTO 100.'

    GOTO 100

  !-----------------------------------------------------------------------
  ! ELSE accept cos_theta_sca and sigma1
  !-----------------------------------------------------------------------     
  ELSE
     IF (chatty) PRINT *,'Accept.'

     cos_theta_sca = ran_cos_theta_sca
     sigma1 = ran_sigma1
        
  ENDIF  

  !-----------------------------------------------------------------------
  ! Safety check
  !-----------------------------------------------------------------------  
  IF (PDF .GT. PDF_max) THEN
     PRINT *,'PDF = ',PDF
     PRINT *,'PDF_max = ',PDF_max
     PRINT *,'Error in rayl_scatter_pol: PDF > PDF_max'
  ENDIF

  !-----------------------------------------------------------------------
  ! 3.) Compute theta (scattered zenith angle), dphi and sigma2
  !     using spherical trigonometry
  !----------------------------------------------------------------------- 
  sin_theta_in = DSQRT(MAX(eps,1.D0-cos_theta_in**2))
  sin_theta_sca = DSQRT(MAX(eps,1.D0-cos_theta_sca**2))

  cos_theta_out = cos_theta_in*cos_theta_sca + sin_theta_in*sin_theta_sca*DCOS(sigma1)
  IF (chatty) PRINT *,'cos_theta_out',cos_theta_out

  sin_theta_out = DSQRT(MAX(eps,1.D0-cos_theta_out**2))

  frac = (cos_theta_sca - cos_theta_in*cos_theta_out) / (sin_theta_in*sin_theta_out)
  IF (frac .GE. (1.D0-eps)) THEN
     dphi = 0.D0
  ELSE IF (frac .LE. (-1.D0+eps)) THEN
     dphi = pi
  ELSE
     dphi = DACOS( frac ) * DSIGN(1.D0,sigma1) 
  ENDIF

  !-----------------------------------------------------------------------
  ! Convert dphi = phi0-phi to dphi = phi-phi0
  !----------------------------------------------------------------------- 
  dphi = -dphi

  !CALL ran(ran3)
  !dphi = 2.D0*pi*ran3

END SUBROUTINE rayl_scatter_pol

!+-----------------------------------------------------------------
! SUBROUTINE mie_scatter_pol_bound
!+-----------------------------------------------------------------
! Sample the cosine of a photon scattering angle and rotation angle
! using rejection sampling,
! and compute scattering zenith angle, azimuth difference,
! and second rotation angle, for Mie scattering.
!
! Victor Trees, September 2022
!+-----------------------------------------------------------------
SUBROUTINE mie_scatter_pol_bound(backward_mode,nmat,PDFerror,                 &
                & new_scatt_data_pol,n_scat_angles,                           &
                & S_in,cos_theta_in,cos_theta_out,dphi,cos_theta_sca)

  !USE knuth_random
  USE new_data_types

  IMPLICIT NONE

  ! Input --------------------------------------------------------
  INTEGER,INTENT(IN) :: nmat !Size of the Stokes vector: 1, 3 or 4 (-)
  INTEGER,INTENT(IN) :: n_scat_angles !Number of single scattering angles (-)
  REAL*8,INTENT(IN)  :: cos_theta_in !Cosine of incident zenith angle (-)
  REAL*8,INTENT(IN)  :: S_in(nmat) !Incident Stokes vector (-)
  LOGICAL            :: backward_mode !Indicates if we run in backward mode

  TYPE(new_scatterer_info_pol),INTENT(IN)  :: new_scatt_data_pol !Object with droplet single scattering info

  ! Output -------------------------------------------------------
  REAL*8,INTENT(OUT) :: cos_theta_out !Cosine of scattered zenith angle (-)
  REAL*8,INTENT(OUT) :: dphi !Azimuth difference between incident and scattered beams (rad)
  REAL*8,INTENT(OUT) :: cos_theta_sca !Cosine of scattering angle

  ! Internal -----------------------------------------------------
  REAL*8 :: scat_angles(2),CDF_bound(2)
  REAL*8 :: sigma1 !polarization rotation angle
  REAL*8 :: PDF,point,frac,PDF_bound,ran_CDF 
  REAL   :: ran1,ran2,ran3
  REAL*8 :: ran_cos_theta_sca,ran_sigma1
  REAL*8 :: sin_theta_in,sin_theta_out,sin_theta_sca
  REAL*8 :: F11,F12,sign

  REAL*8 :: val,low,high,val_mid,val_mid2 
  INTEGER :: goal_found,middle

  REAL*8,PARAMETER :: pi=3.141592653589793D0
  REAL*8,PARAMETER :: EPS=1.0D-16

  LOGICAL             :: chatty,pauses,bruteforce
  LOGICAL             :: PDFerror

  chatty = .false.
  pauses = .false.
  PDFerror = .false.
  bruteforce = .false.

  !-----------------------------------------------------------------------
  ! Sample cos_theta_sca and sigma1
  !-----------------------------------------------------------------------
  IF (chatty) PRINT *,'Sample cos_theta_sca and sigma1'

  !-----------------------------------------------------------------------
  ! Sample random cosine of scattering angle between -1 and 1
  !-----------------------------------------------------------------------

  !Get random CDF value 
202  CALL RANDOM_NUMBER(ran1)  
  !CALL ran(ran1)

  !Convert random number to double precision
  ran_CDF = REAL(ran1,8)

  !Get min and max indices that belong to this random CDF value
  low = 1
  high = n_scat_angles
  goal_found = 0
  val = ran_CDF

  DO WHILE(low <= high .AND. goal_found == 0)
     middle = (low + high)/2

     ! WRITE(*,*)"Now searching element: ", middle

     val_mid = new_scatt_data_pol%CDF_bound(middle)
     val_mid2 = new_scatt_data_pol%CDF_bound(middle+1)

     IF ( (val .GE. val_mid) .AND. (val .LE. val_mid2 ) ) THEN
        goal_found = 1
     ELSE IF ( val .LT. val_mid ) THEN
        high = middle !- 1
     ELSE
        low = middle !+ 1
     END IF
  END DO

  IF (goal_found == 0) THEN
     WRITE(*,*)"Error: ran_CDF not found in CDF_bound."
     STOP
  END IF  

  !Get scattering angles
  scat_angles = new_scatt_data_pol%scat_angles(middle:middle+1)

  !Get cumulative distribution function of bound function for rejection sampling
  CDF_bound = new_scatt_data_pol%CDF_bound(middle:middle+1)

  !sample random cos_theta_sca from CDF_bound
  CALL linpol_dp_single(CDF_bound, scat_angles, 2, ran_CDF, ran_cos_theta_sca)

  IF (chatty) PRINT *,'CDF_bound',CDF_bound
  IF (chatty) PRINT *,'scat_angles',scat_angles
  IF (chatty) PRINT *,'ran_CDF',ran_CDF
  IF (chatty) PRINT *,'ran_cos_theta_sca',ran_cos_theta_sca
  IF (pauses) READ(*,*)

  !-----------------------------------------------------------------------
  ! Sample random rotation angle sigma1 between -pi and pi
  !-----------------------------------------------------------------------
  !CALL ran(ran2)
  CALL RANDOM_NUMBER(ran2)
  ran_sigma1 = 2.D0*pi*ran2 - pi

  IF (chatty) PRINT *,'ran_sigma1',ran_sigma1

  !-----------------------------------------------------------------------
  ! Get Mie scattering matrix for this scattering angle
  !-----------------------------------------------------------------------
  CALL mie_scat_elements(ran_cos_theta_sca,new_scatt_data_pol,n_scat_angles,F11,F12)
  
  IF (nmat .EQ. 1) F12 = 0.D0

  !IF (chatty) PRINT *,'F_sca',F_sca
  IF (chatty) PRINT *,'F11',F11
  IF (chatty) PRINT *,'F12',F12

  !-----------------------------------------------------------------------
  ! Set sign
  !-----------------------------------------------------------------------
  IF (backward_mode) THEN
     sign = 1.D0
  ELSE
     sign = -1.D0      
  ENDIF

  !-----------------------------------------------------------------------
  ! Compute PDF(costheta,sigma1) for Mie scattering
  !-----------------------------------------------------------------------
  IF (nmat .EQ. 1) THEN
     !PDF = F11*S_in(1)
     PDF = F11
  ELSE
     !PDF = F11*S_in(1) + F12*(S_in(2)*DCOS(2.D0*ran_sigma1)-S_in(3)*DSIN(2.D0*ran_sigma1))
     !PDF = F11 + F12*(S_in(2)*DCOS(2.D0*ran_sigma1)-S_in(3)*DSIN(2.D0*ran_sigma1))
     PDF = F11 + F12*(S_in(2)*DCOS(2.D0*ran_sigma1)+sign*S_in(3)*DSIN(2.D0*ran_sigma1))
  ENDIF

  IF (DABS(S_in(1) - 1.D0).GT.eps) THEN
     PRINT *,'S_in=',S_in
     PRINT *,'PDF',PDF
     READ(*,*)
  ENDIF

  IF (nmat .GT. 1) THEN
     !IF (F12*(S_in(2)*DCOS(2.D0*ran_sigma1)-S_in(3)*DSIN(2.D0*ran_sigma1)) .GT. DABS(F12)*DSQRT(2.D0)) THEN
     !IF (F12*(S_in(2)*DCOS(2.D0*ran_sigma1)+sign*S_in(3)*DSIN(2.D0*ran_sigma1)) .GT. DABS(F12)*DSQRT(2.D0)) THEN
     IF (F12*(S_in(2)*DCOS(2.D0*ran_sigma1)+sign*S_in(3)*DSIN(2.D0*ran_sigma1)) .GT. DABS(F12)) THEN
        PRINT *,'S_in(1)=',S_in(1)
        PRINT *,'S_in(2)=',S_in(2)
        PRINT *,'S_in(3)=',S_in(3)
        PRINT *,'PDF',PDF
        READ(*,*)
     ENDIF
  ENDIF  

  !PRINT*,'F_sca=',F_sca

  !-----------------------------------------------------------------------
  ! Generate random point between 0 and max of PDF
  !-----------------------------------------------------------------------
  !CALL ran(ran3)
  CALL RANDOM_NUMBER(ran3)

  !PDF_bound = F11+DABS(F12)*DSQRT(2.D0)
  PDF_bound = F11+DABS(F12)
  
  point = ran3*PDF_bound

  IF (chatty) PRINT *,'point=',point
  IF (pauses) READ(*,*)

  !-----------------------------------------------------------------------
  ! Reject cos_theta_sca and sigma1 if point is larger than PDF
  !-----------------------------------------------------------------------
  IF (point .GT. PDF) THEN

    IF (chatty) PRINT *,'Reject. GOTO 202.'
    IF (pauses) READ(*,*)

    GOTO 202

  !-----------------------------------------------------------------------
  ! ELSE accept cos_theta_sca and sigma1
  !-----------------------------------------------------------------------
  ELSE
     IF (chatty) PRINT *,'Accept.'
     IF (pauses) READ(*,*)

     cos_theta_sca = ran_cos_theta_sca
     sigma1 = ran_sigma1

  ENDIF

  !-----------------------------------------------------------------------
  ! Safety check
  !-----------------------------------------------------------------------
  !IF (PDF .GT. PDF_max) THEN
  IF (PDF .GT. PDF_bound) THEN
     PRINT *,'cos_theta_sca=',cos_theta_sca
     PRINT *,'sigma1=',sigma1
     !PRINT *,'sigma1_try_max=',sigma1_try_max
     PRINT *,'S_in=',S_in
     PRINT *,'F11=',F11
     PRINT *,'F12=',F12
     !PRINT *,'F_sca_max=',F_sca_max
     PRINT *,'PDF = ',PDF
     !PRINT *,'PDF_max = ',PDF_max
     PRINT *,'PDF_bound = ',PDF_bound
     PRINT *,'Error in mie_scatter_pol: PDF > PDF_bound'
     PDFerror = .true.
  ENDIF

  !-----------------------------------------------------------------------
  ! 3.) Compute theta (scattered zenith angle), dphi and sigma2
  !     using spherical trigonometry
  !-----------------------------------------------------------------------
  sin_theta_in = DSQRT(MAX(eps,1.D0-cos_theta_in**2))
  sin_theta_sca = DSQRT(MAX(eps,1.D0-cos_theta_sca**2))

  cos_theta_out = cos_theta_in*cos_theta_sca + sin_theta_in*sin_theta_sca*DCOS(sigma1)
  IF (chatty) PRINT *,'cos_theta_out',cos_theta_out

  sin_theta_out = DSQRT(MAX(eps,1.D0-cos_theta_out**2))

  frac = (cos_theta_sca - cos_theta_in*cos_theta_out) / (sin_theta_in*sin_theta_out)
  IF (frac .GE. (1.D0-eps)) THEN
     dphi = 0.D0
  ELSE IF (frac .LE. (-1.D0+eps)) THEN
     dphi = pi
  ELSE
     dphi = DACOS( frac ) * DSIGN(1.D0,sigma1)
  ENDIF

  !-----------------------------------------------------------------------
  ! Convert dphi = phi0-phi to dphi = phi-phi0
  !-----------------------------------------------------------------------
  dphi = -dphi

  !CALL ran(ran3)
  !dphi = 2.D0*pi*ran3

END SUBROUTINE mie_scatter_pol_bound


!+-----------------------------------------------------------------
! SUBROUTINE stokes_rot_mat
!+-----------------------------------------------------------------
! Returns the Stokes transformation matrices 
!
! Victor Trees
!+-----------------------------------------------------------------
SUBROUTINE stokes_rot_mat(L1,L2,cos_theta_sca,uo,u1,phio,phi1,nmat)
     
  IMPLICIT NONE

  ! Input --------------------------------------------------------
  INTEGER,INTENT(IN) :: nmat          !Size of the Stokes vector: 1, 3 or 4 (-)
  REAL*8,INTENT(IN)  :: cos_theta_sca !Cosine of single scattering angle (-)
  REAL*8,INTENT(IN)  :: uo,u1         !Cosines of incident (scattered) zenith angles (-)
  REAL*8,INTENT(IN)  :: phio          !Azimuth angle of inicident photon (rad) 
  REAL*8,INTENT(IN)  :: phi1          !Azimuth angle of scattered photon (rad)

  ! Output -------------------------------------------------------
  REAL*8,DIMENSION(nmat,nmat),INTENT(OUT)    :: L1,L2 !Stokes transformation matrices (-)

  ! Internal -----------------------------------------------------
  INTEGER                              :: casetype 
  REAL*8                               :: cosi1,cosi2,qv(3),q,sini1,sini2
  REAL*8                               :: wo(3),w1(3),wophi(3),w1phi(3)
  REAL*8                               :: sino,sin1,fac,dphi
  REAL*8                               :: sin_theta_sca
  REAL*8,PARAMETER                     :: eps=1.0D-16 
  REAL*8,PARAMETER                     :: DegtoRad=0.017453292519943295D0
  REAL*8,PARAMETER                     :: pi=3.141592653589793D0
  LOGICAL                              :: chatty, pauses
  
  chatty = .false.
  pauses = .false.

  !-----------------------------------------------------------------------
  ! Initialize rotation matrices
  !-----------------------------------------------------------------------   
  L1=0.0D+0
  L2=0.0D+0

  !-----------------------------------------------------------------------
  ! Set diagonal corners of rotation matrices
  !-----------------------------------------------------------------------   
  L1(1,1)=1.0D+0
  L2(1,1)=1.0D+0

  IF (nmat.LT.3) THEN
     GOTO 911
  ENDIF

  IF (nmat.GT.3) THEN
     L1(4,4)=1.0D+0
     L2(4,4)=1.0D+0
  ENDIF

  !-----------------------------------------------------------------------
  ! If meridian planes coincide, no rotation is needed.
  !-----------------------------------------------------------------------
  dphi = DABS(phio-phi1)
  
  !IF (ABS(phi1-phio).LT.eps) THEN

  IF ((DMOD(dphi,pi) .LT. eps).OR.(DABS(DMOD(dphi,pi)-pi).LT.eps)) THEN
  
    casetype = 1      

    L1(2,2)=1.0D+0
    L1(3,3)=1.0D+0

    L2(2,2)=1.0D+0
    L2(3,3)=1.0D+0

    !PRINT *,'Both scattered and incident beam are vertical!'

    IF (chatty) WRITE(*,*) 'L1 = ',L1
    IF (chatty) WRITE(*,*) 'L2 = ',L2
    IF (pauses) READ(*,*)

    GOTO 911

  ENDIF

  !-----------------------------------------------------------------------
  ! Skip if both scattered and incident beams are vertical
  ! (Case 4 of Hovenier et al., 2004).
  !-----------------------------------------------------------------------
  IF ((DABS(uo).GT.1.0D+0-eps).AND.(DABS(u1).GT.1.0D+0-eps)) THEN

    casetype = 2

    L1(2,2)=1.0D+0
    L1(3,3)=1.0D+0

    L2(2,2)=1.0D+0
    L2(3,3)=1.0D+0

    GOTO 911

  ENDIF

  !-----------------------------------------------------------------------
  ! Compute wo (incident vector) and w1 (scattered vector)
  !-----------------------------------------------------------------------
  fac = 1.0D+0-uo**2
  IF (fac .LT. eps) THEN
     sino = 0.D0
  ELSE
     sino=DSQRT(fac)
  ENDIF

  fac = 1.0D+0-u1**2
  IF (fac .LT. eps) THEN
     sin1 = 0.D0
  ELSE
     sin1=DSQRT(fac)
  ENDIF

  wo=[sino*DCOS(phio),sino*DSIN(phio),uo] !x,y,z in pol scat plane
  w1=[sin1*DCOS(phi1),sin1*DSIN(phi1),u1] !x,y,z in pol scat plane

  !-----------------------------------------------------------------------
  ! IF special cases (incident or scattered beam close to vertical)
  !-----------------------------------------------------------------------
  IF ((DABS(uo).GT.1.0D+0-eps).OR.(DABS(u1).GT.1.0D+0-eps)) THEN

    !-----------------------------------------------------------------------
    ! Define azimuth vectors on horizontal plane
    !-----------------------------------------------------------------------          
    wophi=[DCOS(phio),DSIN(phio),0.0D+0]
    w1phi=[DCOS(phi1),DSIN(phi1),0.0D+0]

    !above two lines are similar to:
    !wophi=[sino*DCOS(phio),sino*DSIN(phio),0.0D+0]
    !w1phi=[sin1*DCOS(phi1),sin1*DSIN(phi1),0.0D+0] 
    !if you divied by sino and sin1

    !-----------------------------------------------------------------------
    ! IF scattered beam close to vertical (up or down)
    !-----------------------------------------------------------------------
    !if (abs(uo).gt.abs(u1)) then
    IF (DABS(u1).GT.1.0D+0-eps) THEN                          

      casetype = 3

      cosi1=DSIGN(1.0D+0,u1)                                   
      cosi2=-cosi1*(w1phi(1)*wophi(1)+w1phi(2)*wophi(2))

      qv=[wo(2),-wo(1),0.0D+0]
      q=qv(1)*w1phi(1)+qv(2)*w1phi(2)

    !-----------------------------------------------------------------------
    ! IF incident beam close to vertical (up or down)
    !-----------------------------------------------------------------------
    ELSE IF (DABS(uo).GT.1.0D+0-eps) THEN                    

      casetype = 4      

      cosi2=DSIGN(1.0D+0,uo)                                 
      cosi1=-cosi2*(w1phi(1)*wophi(1)+w1phi(2)*wophi(2))

      qv=[wophi(2),-wophi(1),0.0D+0]
      q=qv(1)*w1(1)+qv(2)*w1(2)

    ENDIF

  !-----------------------------------------------------------------------
  ! ELSE there is NO special case
  !-----------------------------------------------------------------------    
  ELSE

    casetype = 5

    cosi1=(u1-uo*cos_theta_sca)/DSQRT(MAX(eps,(1.0D+0-uo**2)*(1.0D+0-cos_theta_sca**2))) !<<< good
    cosi2=(uo-u1*cos_theta_sca)/DSQRT(MAX(eps,(1.0D+0-u1**2)*(1.0D+0-cos_theta_sca**2))) !<<< good

    qv=[wo(2),-wo(1),0.0D+0]
    q=qv(1)*w1(1)+qv(2)*w1(2)

  ENDIF

  !-----------------------------------------------------------------------
  ! Compute sini1 and sini2
  !-----------------------------------------------------------------------   
  fac = 1.0D+0-cosi1**2
  IF (fac .LT. eps) THEN
     sini1 = 0.D0
  ELSE
     sini1=DSQRT(fac)*DSIGN(1.D0,q)
  ENDIF  

  fac = 1.0D+0-cosi2**2
  IF (fac .LT. eps) THEN
     sini2 = 0.D0
  ELSE
     sini2=DSQRT(fac)*DSIGN(1.D0,q)
  ENDIF

  !-----------------------------------------------------------------------
  ! Compute elements rotation matrices
  !----------------------------------------------------------------------- 
  L1(2,2)=cosi1**2-sini1**2
  L1(3,3)=L1(2,2)
  L1(2,3)=-2.0D+0*sini1*cosi1
  L1(3,2)=-L1(2,3)

  L2(2,2)=cosi2**2-sini2**2
  L2(3,3)=L2(2,2)
  L2(2,3)=-2.0D+0*sini2*cosi2
  L2(3,2)=-L2(2,3)

  IF (chatty) THEN
     IF ((DABS(L1(2,2)).GT.1.0D0+1.0D-10) .OR. (DABS(L1(2,3)).GT.1.0D0+1.0D-10) .OR. &
          & (DABS(L2(2,2)).GT.1.0D0+1.0D-10) .OR. (DABS(L2(2,3)).GT.1.0D0+1.0D-10) ) THEN
        PRINT *,'L1(2,2)=',L1(2,2)
        PRINT *,'L1(2,3)=',L1(2,3)
        PRINT *,'L2(2,2)=',L2(2,2)
        PRINT *,'L2(2,3)=',L2(2,3)
        PRINT *,'casetype=',casetype
        PRINT *,'sino=',sino,', sin1=',sin1
        PRINT *,'cos_theta_sca=',cos_theta_sca
        !PRINT *,'cos_theta_sca_2=',cos_theta_sca_2
        PRINT *,'sin_theta_sca=',sin_theta_sca
        PRINT *,'cosi1=',cosi1,', sini1=',sini1
        PRINT *,'cosi2=',cosi2,', sini2=',sini2
        PRINT *,'wo=',wo
        PRINT *,'w1=',w1
        PRINT *,'wophi=',wophi
        PRINT *,'w1phi=',w1phi
        PRINT *,'q=',q
        PRINT *,'qv=',qv
        PRINT *,'phio=',phio,', phi1=',phi1
        PRINT *,'uo=',uo,', u1=',u1
        READ(*,*)
     ENDIF
  ENDIF

911 RETURN

  END SUBROUTINE stokes_rot_mat

!+-----------------------------------------------------------------
! SUBROUTINE matmul_stokes_rot
!+-----------------------------------------------------------------
! Left and right multiply Stokes matrix by polarization rotation
! matrices L1 and L2.
!
! Victor Trees
!+-----------------------------------------------------------------
SUBROUTINE matmul_stokes_rot(nmat,A,L1,L2,Arot)

  IMPLICIT NONE

  ! Input --------------------------------------------------------
  INTEGER,INTENT(IN)  :: nmat !Size of the Stokes vector: 1, 3 or 4 (-)
  REAL*8,INTENT(IN)   :: A(nmat,nmat) !Single scattering matrix matrix (-)
  REAL*8,INTENT(IN)   :: L1(nmat,nmat),L2(nmat,nmat) !Rotation matrices (-)

  ! Output -------------------------------------------------------
  REAL*8,INTENT(OUT)   :: Arot(nmat,nmat) !Phase matrix (-)

  ! Internal -----------------------------------------------------
  REAL*8 :: C1,S1,C2,S2

!-----------------------------------------------------------------------
! Left and right multiply matrix A with the rotation matrices
! analogous to Eq. (88) of Hovenier and van der Mee (1983).
! Note that Hovenier and van der Mee chose to substitute their
! scattering matrix elements (a1, b1, b2, -b2,...) already, but we
! keep the entries as is (i.e. no explicit minus sign for A(4,3)):
!-----------------------------------------------------------------------
  Arot(1,1) = A(1,1)

  IF (nmat.GT.1) THEN

     C1 = L1(2,2)
     S1 = L1(3,2)
     C2 = L2(2,2)
     S2 = L2(3,2)

     Arot(1,2) = A(1,2)*C1
     Arot(1,3) = -A(1,2)*S1

     Arot(2,1) = C2*A(2,1)
     Arot(2,2) = C1*C2*A(2,2) - S1*S2*A(3,3)
     Arot(2,3) = -C2*S1*A(2,2) - C1*S2*A(3,3)

     Arot(3,1) = S2*A(2,1)
     Arot(3,2) = C1*S2*A(2,2) + C2*S1*A(3,3)
     Arot(3,3) = C1*C2*A(3,3) - S1*S2*A(2,2)

     IF (nmat.EQ.4) THEN

        Arot(4,1) = 0.D0
        Arot(4,2) = S1*A(4,3)
        Arot(4,3) = C1*A(4,3)
        Arot(4,4) = A(4,4)

        Arot(1,4) = 0.D0
        Arot(2,4) = -S2*A(3,4)
        Arot(3,4) = C2*A(3,4) 

     ENDIF

  ENDIF 

  RETURN
END SUBROUTINE matmul_stokes_rot


END MODULE mc_routines_general


