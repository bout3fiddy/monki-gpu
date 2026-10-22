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
! MODULE mc_routines_cartesian
!
! Module that contains subroutines for the mc.f90 program of MONKI.
! Subroutines that only apply to cartesian version of the code. 
!
! Author: Victor Trees
!------------------------------------------------------------------
MODULE mc_routines_cartesian

  IMPLICIT NONE

CONTAINS

!+-----------------------------------------------------------------
! SUBROUTINE calc_radiances_pol
!+-----------------------------------------------------------------
! Compute reflectance contribution of this scattering event by
! by propagating the scattered light towards the detector.
! Version without separate scattering contributions 
! or air mass factors.
!
! Victor Trees
!+-----------------------------------------------------------------
SUBROUTINE calc_radiances_pol(nx, ny, nz, nmat, n_angles, &
     &     a, b, c, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
     &     TAU_gas,lowest_k_of_only_gas, &
     &     blocked_sensor_path, precompute_detector_paths, &
     &     va, vb, vc,k_detector,&
     &     S_prob,Stokes,k_register,backward_mode,i_inject,j_inject)

  IMPLICIT NONE

  ! Input -------------------------------------------------------------------------------------------
  INTEGER,INTENT(IN) :: nx,ny,nz              !Number of grid cells in x, y and z directions (-)
  INTEGER,INTENT(IN) :: nmat                  !Size of the Stokes vector: 1, 3 or 4 (-)
  INTEGER,INTENT(IN) :: n_angles              !Number of viewing zenith angles (-)
  INTEGER,INTENT(IN) :: k_detector            !Layer index of the detector location (-)
  INTEGER,INTENT(IN) :: k_register            !Layer index of the save grid (-)
  INTEGER,INTENT(IN) :: lowest_k_of_only_gas  !Lowest layer index above which only gas occurs (-)
  INTEGER,INTENT(IN) :: i_inject,j_inject     !Column indices of photon injection in backward mode (-)
  INTEGER,INTENT(IN) :: blocked_sensor_path(nx,ny,nz,n_angles) !Indicates if path to detector is blocked by aerosol (-)

  REAL*8,INTENT(IN) :: BETA(NX, NY, NZ)        !Extinction coefficient (1/km)
  REAL*8,INTENT(IN) :: xunit, yunit, ZUNIT(NZ) !Dimensions of the grid cells in the x, y and z diretion (km)
  REAL*8,INTENT(IN) :: VA(N_ANGLES), VB(N_ANGLES), VC(N_ANGLES) !x, y and z direction towards detector (-)
  REAL*8,INTENT(IN) :: S_PROB(N_ANGLES,nmat)   !Reflected Stokes vector at location of event (-)
  REAL*8,INTENT(IN) :: TAU_gas(NZ,N_ANGLES)    !Optical depth of path through gas towards detector (-)

  LOGICAL,INTENT(IN) :: backward_mode          !Indicates if the code runs in backward mode
  LOGICAL,INTENT(IN) :: precompute_detector_paths !Indicate if detector paths are precomputed

  ! In- and output ---------------------------------------------------------------------------------
  INTEGER,INTENT(INOUT) :: i,j,k               !Indices of the x, y and z grid cells

  REAL*8,INTENT(INOUT) :: a,b,c                !x, y and z directions of the photon before scattering (-)
  REAL*8,INTENT(INOUT) :: x,y,z                !Straight distances from photon to next x, y and z boundaries (km)
  REAL*8,INTENT(INOUT) :: STOKES(nmat,NX, NY, N_ANGLES) !Accumulated signal to the detector (-)

  ! Internal --------------------------------------------------------------------------------------
  INTEGER :: DELTA_I,DELTA_J,I_NEW,J_NEW,direction
  INTEGER :: i_bak,j_bak,k_bak,ia
  INTEGER :: i_sign,j_sign,k_sign
  INTEGER :: K_OLD

  REAL*8 :: tau_v,dist, beta_value
  REAL*8 :: a_bak,b_bak,c_bak
  REAL*8 :: a_abs,b_abs,c_abs,a_inv,b_inv,c_inv
  REAL*8 :: x_bak,y_bak,z_bak
  REAL*8 :: HEIGHT,DELTA_X,DELTA_Y,x_old,y_old,z_old
  REAL*8 :: CONTRIB_STOKES(nmat),CONTRIB_STOKES_NORM(nmat)

  REAL*8,PARAMETER :: eps = 1.0D-16

  LOGICAL :: chatty,pauses,advance_test

  chatty = .false.      !Print comments to screen
  pauses = .false.      !Pause script regularly
  advance_test = .false. !Double check if projection is okay

  !+-----------------------------------------------------------------
  ! Record actual photon position and direction.
  !+-----------------------------------------------------------------
  A_BAK = A
  B_BAK = B
  C_BAK = C
  I_BAK = I
  J_BAK = J
  K_BAK = K
  X_BAK = X
  Y_BAK = Y
  Z_BAK = Z

  !+-----------------------------------------------------------------
  ! Loop over each viewing direction.  Reset the photon position and 
  ! direction.  Zero TAU_V, the optical thickness remaining in the 
  ! model in the given direction.
  !+-----------------------------------------------------------------
  DO IA=1,N_ANGLES
     IF (CHATTY) PRINT*,'IA=',IA

     A = VA(IA)
     B = VB(IA)
     C = VC(IA)
     I = I_BAK
     J = J_BAK
     K = K_BAK
     X = X_BAK
     Y = Y_BAK
     Z = Z_BAK

     !+-----------------------------------------------------------------
     ! If the scattered photon travels towards the sensor, then:
     !+-----------------------------------------------------------------
     !IF (((k.LE.k_detector).AND.(vc(ia).GT.0.D0)).OR.((k.GE.k_detector).AND.(vc(ia).LT.0.D0))) THEN
     !   IF (CHATTY) PRINT*,'Scattered photon travels towards sensor'

        IF (S_PROB(IA,1) .GT. EPS) THEN

           IF (CHATTY) PRINT*,'S_PROB(IA,1)=',S_PROB(IA,1)
           IF (CHATTY) PRINT*,'(S_PROB(IA,1) .GT. EPS)'

           !+-----------------------------------------------------------------
           ! If the sign of direction changed, simply reverse x, y and/or z
           !+-----------------------------------------------------------------
           IF (A*A_BAK .LT. 0.D0) X = MAX(0.D0, XUNIT - X)
           IF (B*B_BAK .LT. 0.D0) Y = MAX(0.D0, YUNIT - Y)
           IF (C*C_BAK .LT. 0.D0) Z = MAX(0.D0, ZUNIT(K) - Z)

           IF (CHATTY) PRINT*,'X=',X,'Y=',Y,'Z=',Z

           !+-----------------------------------------------------------------
           ! Zero TAU_V and path
           !+-----------------------------------------------------------------
           TAU_V      = 0.D0

           !+-----------------------------------------------------------------
           ! Get signs, abs values and inverses
           !+-----------------------------------------------------------------
           I_SIGN     = INT(SIGN(1.D0, A))
           J_SIGN     = INT(SIGN(1.D0, B))
           K_SIGN     = INT(SIGN(1.D0, C))
           A_ABS      = ABS(A)
           B_ABS      = ABS(B)
           C_ABS      = ABS(C)
           A_INV      = 1.D0 / MAX(A_ABS, EPS)
           B_INV      = 1.D0 / MAX(B_ABS, EPS)
           C_INV      = 1.D0 / MAX(C_ABS, EPS)

           !+-----------------------------------------------------------------
           ! Store photon location, with signs in direction towards sensor
           !+-----------------------------------------------------------------           
           x_old = x
           y_old = y
           z_old = z

           !+-----------------------------------------------------------------
           ! Compute TAU_V from event to sensor
           !+----------------------------------------------------------------- 
           IF (CHATTY) PRINT*,'K=',K,'NZ=',NZ,'K_detector=',K_detector

           !As long as the photon is still in the scene
           !DO WHILE ((K .GE. 1 .AND. K .LE. NZ).AND.(k .NE. k_detector))
           DO WHILE (K .GE. 1 .AND. K .LT. NZ)

              !!and if the photon is above the clouds             
              !IF (k.GE.lowest_k_of_only_gas) THEN

              !and if we do not save the path (only then we can do this)
              !IF (.NOT. save_path) THEN

              IF (precompute_detector_paths) THEN

                 !and if the path to the sensor is not blocked by aerosol            
                 IF (blocked_sensor_path(I,J,K,IA) .EQ. 0) THEN
                 !IF (.false.) THEN

                    !and if the photon is at a z-boundary
                    IF (DABS(Z-ZUNIT(K)).LE.EPS) THEN

                       !then we can just use TAU_gas
                       TAU_V = TAU_V + TAU_gas(K,IA)

                       k = k_detector !this line is not really needed

                       GOTO 5009

                    ENDIF
                 ENDIF
              ENDIF

              !Otherwise we compute the path optical thickness per cell
              BETA_VALUE = BETA(I,J,K)

              k_old = k

              CALL ADVANCE(NX, NY, NZ, I, J, K, X, Y, Z, &
                   &              I_SIGN, J_SIGN, K_SIGN, A_ABS, B_ABS, C_ABS, &
                   &              A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT, DIST)

              TAU_V = TAU_V + DIST * BETA_VALUE

5009       ENDDO

           !-----------------------------------------------------------------------------------
           ! Now see where this event would be registered!
           !-----------------------------------------------------------------------------------
           IF (CHATTY) PRINT*,'See where this event would be registered.'

           !In backward monte carlo mode, we register at the injection location
           IF (backward_mode) THEN
              i = i_inject
              j = j_inject
              GOTO 5010
           ENDIF

           !Check if the detector is in the same direction as the register grid, as seen from the event
           IF (k_register .EQ. k_detector) THEN
              direction = 1
           ELSE
              direction = -1
           ENDIF

           !Double calculation of the mapping for debugging
           IF (advance_test) THEN

              I = I_BAK
              J = J_BAK
              K = K_BAK

              IF (direction .EQ. -1) THEN
                 X = XUNIT-X_OLD
                 Y = YUNIT-Y_OLD
                 Z = ZUNIT(K)-Z_OLD
              ELSE
                 X = X_OLD
                 Y = Y_OLD
                 Z = Z_OLD
              ENDIF

              !DO WHILE ((K .GE. 1 .AND. K .LE. NZ).AND.(k .NE. k_detector))
              DO WHILE (K .GE. 1 .AND. K .LT. NZ) !k = nz is vacuum layer

                 CALL ADVANCE(NX, NY, NZ, I, J, K, X, Y, Z, &
                & direction*I_SIGN, direction*J_SIGN, direction*K_SIGN, A_ABS, B_ABS, C_ABS, &
                &              A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT, DIST)

              ENDDO

           ENDIF

           !Map the signal onto the nx*ny output register grid
           IF (K_BAK .GT. 1) THEN
              HEIGHT = SUM(ZUNIT(1:K_BAK-1))
           ELSE
              HEIGHT = 0.D0
           ENDIF

           IF (K_SIGN .GT. 0) THEN
              HEIGHT = HEIGHT+ZUNIT(K_BAK)-Z_OLD
           ELSE
              HEIGHT = HEIGHT+Z_OLD
           ENDIF

           IF (k_register .EQ. NZ) THEN
              HEIGHT = SUM(ZUNIT(1:NZ-1)) - HEIGHT
           ENDIF

           DELTA_X = A_ABS * HEIGHT / C_ABS
           DELTA_Y = B_ABS * HEIGHT / C_ABS

           IF (direction .EQ. -1) THEN
              DELTA_I = INT( ( X_OLD + DELTA_X ) / XUNIT )
              DELTA_J = INT( ( Y_OLD + DELTA_Y ) / YUNIT )
           ELSE
              DELTA_I = INT( ( (XUNIT-X_OLD) + DELTA_X ) / XUNIT )
              DELTA_J = INT( ( (YUNIT-Y_OLD) + DELTA_Y ) / YUNIT )
           ENDIF

           I_NEW = MODULO(I_BAK + direction*I_SIGN*DELTA_I + NX - 1, NX) + 1
           J_NEW = MODULO(J_BAK + direction*J_SIGN*DELTA_J + NY - 1, NY) + 1

           !Check if the calculation was correction
           IF (advance_test) THEN

           IF ((I_NEW .NE. I) .OR. (J_NEW .NE. J)) THEN
              IF ((I_NEW .NE. I) .AND. (MOD(ABS(X_OLD),XUNIT) .LT. eps)) THEN
                 PRINT *,'Warning: (I_NEW .NE. I): I_NEW=',I_NEW,'I=',I
                 PRINT *,'MOD(ABS(X_OLD),XUNIT)=',MOD(ABS(X_OLD),XUNIT)
              ELSEIF ((J_NEW .NE. J) .AND. (MOD(ABS(Y_OLD),YUNIT) .LT. eps)) THEN
                 PRINT *,'Warning: (J_NEW .NE. J): J_NEW=',J_NEW,'J=',J
                 PRINT *,'MOD(ABS(Y_OLD),YUNIT)=',MOD(ABS(Y_OLD),YUNIT)
              ELSE
                 PRINT *,'ERROR!'
                 PRINT *,'NX=',NX
                 PRINT *,'-1*I_SIGN*DELTA_I=',-1*I_SIGN*DELTA_I
                 PRINT *,'I_BAK + -1*I_SIGN*DELTA_I + NX - 1=',I_BAK - I_SIGN*DELTA_I + NX - 1
                 PRINT *,'MOD(I_BAK + -1*I_SIGN*DELTA_I + NX - 1, NX) + 1',MOD(I_BAK - I_SIGN*DELTA_I + NX - 1, NX) + 1
                 PRINT *,'MOD(-1,20)=',MODULO(-1,20)
                 PRINT *,'.........'
                 PRINT *,'X=',X,'Y=',Y,'Z=',Z
                 PRINT *,'HEIGHT=',HEIGHT
                 PRINT *,'A_ABS=',A_ABS
                 PRINT *,'B_ABS=',B_ABS
                 PRINT *,'C_ABS=',C_ABS
                 PRINT *,'XUNIT=',XUNIT
                 PRINT *,'YUNIT=',YUNIT
                 PRINT *,'X_OLD=',X_OLD
                 PRINT *,'Y_OLD=',Y_OLD
                 PRINT *,'Z_OLD=',Z_OLD
                 PRINT *,'(XUNIT-X_OLD)=',(XUNIT-X_OLD)
                 PRINT *,'(YUNIT-Y_OLD)=',(YUNIT-Y_OLD)
                 PRINT *,'X_BAK=',X_BAK
                 PRINT *,'Y_BAK=',Y_BAK
                 PRINT *,'DELTA_X=',DELTA_X
                 PRINT *,'DELTA_Y=',DELTA_Y
                 PRINT *,'I_BAK=',I_BAK
                 PRINT *,'J_BAK=',J_BAK
                 PRINT *,'K_BAK=',K_BAK
                 PRINT *,'DELTA_I=',DELTA_I
                 PRINT *,'DELTA_J=',DELTA_J
                 PRINT *,'I_SIGN=',I_SIGN
                 PRINT *,'J_SIGN=',J_SIGN
                 PRINT *,'K_SIGN=',K_SIGN
                 PRINT *,'direction=',direction
                 PRINT *,'SUM(ZUNIT(1:NZ-1))=',SUM(ZUNIT(1:NZ-1))
                 PRINT *,'SUM(ZUNIT(1:K_BAK-1))=',SUM(ZUNIT(1:K_BAK-1))
                 PRINT *,'k_register=',k_register
                 PRINT *,'k_detector=',k_detector
                 PRINT *,'NZ=',NZ
                 PRINT *,'I_NEW=',I_NEW,', I=',I
                 PRINT *,'J_NEW=',J_NEW,', J=',J
                 READ(*,*)
              ENDIF
           ENDIF
           ENDIF

           !Update I and J
           I = I_NEW
           J = J_NEW

           !+-----------------------------------------------------------------
           ! Compute contribution
           !+-----------------------------------------------------------------
5010       IF (CHATTY) PRINT*,'Compute contribution'

           CONTRIB_STOKES(:) = S_PROB(IA,:) * EXP(-TAU_V)

           !+-----------------------------------------------------------------
           ! Compute radiances and Stokes
           !+-----------------------------------------------------------------              

           IF (chatty) PRINT*,'STOKES(:,I,J,IA)=',STOKES(:,I,J,IA),'CONTRIB_STOKES(:)',CONTRIB_STOKES(:)

           STOKES(:,I,J,IA) = STOKES(:,I,J,IA) + CONTRIB_STOKES(:)

           IF (chatty) PRINT*,'STOKES(:,I,J,IA)=',STOKES(:,I,J,IA),'CONTRIB_STOKES(:)',CONTRIB_STOKES(:)
           IF (pauses) READ(*,*)

           !ENDIF
        ENDIF
     !ENDIF

  ENDDO

  !+-----------------------------------------------------------------
  ! Restore actual photon position and direction.
  !+-----------------------------------------------------------------
  A = A_BAK
  B = B_BAK
  C = C_BAK
  I = I_BAK
  J = J_BAK
  K = K_BAK
  X = X_BAK
  Y = Y_BAK
  Z = Z_BAK

  RETURN

  END SUBROUTINE CALC_RADIANCES_POL


!+-----------------------------------------------------------------
! SUBROUTINE calc_radiances_and_scats
!+-----------------------------------------------------------------
! Compute reflectance contribution of this scattering event by
! by propagating the scattered light towards the detector.
! Version that saves single, multiple and 3D contributions.
!
! Victor Trees
!+-----------------------------------------------------------------
SUBROUTINE calc_radiances_and_scats(nx, ny, nz, nmat, n_angles, &
     &     nx_1, ny_1, nz_1, nmat_1, n_angles_1, &
     &     a, b, c, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
     &     TAU_gas,lowest_k_of_only_gas, &
     &     blocked_sensor_path, precompute_detector_paths,&
     &     va, vb, vc,k_detector,&
     &     S_prob,Stokes,nscat,&
     &     CONTRIB_STOKES_3D,&
     &     SINGSCAT_CONTRIB_STOKES_3D,&
     &     MULTSCAT_CONTRIB_STOKES_3D,&
     &     Stokes_sing,Stokes_mult, &
     &     k_register,backward_mode,i_inject,j_inject, &
     &     i_store_singscat,j_store_singscat,k_store_singscat)

  IMPLICIT NONE

  ! Input -------------------------------------------------------------------------------------------
  INTEGER,INTENT(IN) :: nx,ny,nz              !Number of grid cells in x, y and z directions (-)
  INTEGER,INTENT(IN) :: nmat                  !Size of the Stokes vector: 1, 3 or 4 (-)
  INTEGER,INTENT(IN) :: n_angles              !Number of viewing zenith angles (-)
  INTEGER,INTENT(IN) :: nx_1,ny_1,nz_1        !Similar, but only relevant if save_scats=.true.
  INTEGER,INTENT(IN) :: nmat_1                !Similar, but only relevant if save_scats=.true.
  INTEGER,INTENT(IN) :: n_angles_1            !Similar, but only relevant if save_scats=.true.
  INTEGER,INTENT(IN) :: k_detector            !Layer index of the detector location (-)
  INTEGER,INTENT(IN) :: k_register            !Layer index of the save grid (-)
  INTEGER,INTENT(IN) :: lowest_k_of_only_gas  !Lowest layer index above which only gas occurs (-)
  INTEGER,INTENT(IN) :: nscat                 !Scattering order counter (-)
  INTEGER,INTENT(IN) :: i_inject,j_inject     !Column indices of photon injection in backward mode (-)

  INTEGER,INTENT(IN) :: blocked_sensor_path(nx,ny,nz,n_angles) !Indicates if path to detector is blocked by aerosol (-)

  REAL*8,INTENT(IN) :: BETA(NX, NY, NZ)        !Extinction coefficient (1/km)
  REAL*8,INTENT(IN) :: xunit, yunit, ZUNIT(NZ) !Dimensions of the grid cells in the x, y and z diretion (km)
  REAL*8,INTENT(IN) :: VA(N_ANGLES), VB(N_ANGLES), VC(N_ANGLES) !x, y and z direction towards detector (-)
  REAL*8,INTENT(IN) :: S_PROB(N_ANGLES,nmat)   !Reflected Stokes vector at location of event (-)
  REAL*8,INTENT(IN) :: TAU_gas(NZ,N_ANGLES)    !Optical depth of path through gas towards detector (-)

  LOGICAL,INTENT(IN) :: backward_mode          !Indicates if the code runs in backward mode
  LOGICAL,INTENT(IN) :: precompute_detector_paths !Indicate if detector paths are precomputed

  ! In- and output ---------------------------------------------------------------------------------
  INTEGER,INTENT(INOUT) :: i,j,k               !Indices of the x, y and z grid cells
  INTEGER,INTENT(INOUT) :: i_store_singscat    !x index of single scattering of photon
  INTEGER,INTENT(INOUT) :: j_store_singscat    !y index of single scattering of photon
  INTEGER,INTENT(INOUT) :: k_store_singscat    !z index of single scattering of photon

  REAL*8,INTENT(INOUT) :: a,b,c                !x, y and z directions of the photon before scattering (-)
  REAL*8,INTENT(INOUT) :: x,y,z                !Straight distances from photon to next x, y and z boundaries (km)
  REAL*8,INTENT(INOUT) :: STOKES(nmat,NX, NY, N_ANGLES) !Accumulated signal to the detector (-)
  REAL*8,INTENT(INOUT) :: STOKES_SING(nmat_1,NX_1, NY_1, N_ANGLES_1) !Signal of single scattering only (-)
  REAL*8,INTENT(INOUT) :: STOKES_MULT(nmat_1,NX_1, NY_1, N_ANGLES_1) !Signal of multiple scattering only (-)
  REAL*8,INTENT(INOUT) :: CONTRIB_STOKES_3D(nmat_1,NX_1,NY_1,0:NZ_1,N_ANGLES_1) !Layer conribution to signal (-)
  REAL*8,INTENT(INOUT) :: SINGSCAT_CONTRIB_STOKES_3D(nmat_1,NX_1,NY_1,0:NZ_1,N_ANGLES_1) !Layer contrib to sing scat signal (-)
  REAL*8,INTENT(INOUT) :: MULTSCAT_CONTRIB_STOKES_3D(nmat_1,NX_1,NY_1,0:NZ_1,N_ANGLES_1) !Layer contrib to mult scat signal (-)

  ! Internal --------------------------------------------------------------------------------------
  INTEGER :: DELTA_I,DELTA_J,I_NEW,J_NEW,direction
  INTEGER :: i_bak,j_bak,k_bak,ia
  INTEGER :: i_sign,j_sign,k_sign
  INTEGER :: I_OLD,J_OLD,K_OLD,k_store,i_store,j_store

  REAL*8 :: tau_v,dist, beta_value
  REAL*8 :: a_bak,b_bak,c_bak
  REAL*8 :: a_abs,b_abs,c_abs,a_inv,b_inv,c_inv
  REAL*8 :: x_bak,y_bak,z_bak
  REAL*8 :: HEIGHT,DELTA_X,DELTA_Y,x_old,y_old,z_old
  REAL*8 :: CONTRIB_STOKES(nmat),CONTRIB_STOKES_NORM(nmat)

  REAL*8,PARAMETER :: eps = 1.0D-16

  LOGICAL :: chatty,pauses,advance_test

  chatty = .false.      !Print comments to screen
  pauses = .false.      !Pause script regularly
  advance_test = .false. !Double check if projection is okay

  !+-----------------------------------------------------------------
  ! Record actual photon position and direction.
  !+-----------------------------------------------------------------
  A_BAK = A
  B_BAK = B
  C_BAK = C
  I_BAK = I
  J_BAK = J
  K_BAK = K
  X_BAK = X
  Y_BAK = Y
  Z_BAK = Z

  !+-----------------------------------------------------------------
  !Check if this was a surface reflection: then k_store = 0
  !+-----------------------------------------------------------------
  IF ((k_bak .EQ. 1) .AND. (z_bak .LT. eps)) THEN
     k_store = 0
  ELSE
     k_store = k_bak
  ENDIF

  !+-----------------------------------------------------------------
  !Check if this was first order scattering in backward mode,
  !if yes, save location, because in the backward mode the sensor
  !is directly looking at the location of 1st order scattering
  !(instead of the last order in forward mode)
  !+-----------------------------------------------------------------
  IF (backward_mode) THEN
     IF (nscat .EQ. 1) THEN
        i_store = i_bak
        j_store = j_bak
        i_store_singscat = i_store
        j_store_singscat = j_store
        k_store_singscat = k_store
     ELSE
        i_store = i_store_singscat
        j_store = j_store_singscat
        k_store = k_store_singscat
     ENDIF
  ELSE
     i_store = i_bak
     j_store = j_bak
  ENDIF

  !+-----------------------------------------------------------------
  ! Loop over each viewing direction.  Reset the photon position and
  ! direction.  Zero TAU_V, the optical thickness remaining in the
  ! model in the given direction.
  !+-----------------------------------------------------------------
  DO IA=1,N_ANGLES
     IF (CHATTY) PRINT*,'IA=',IA

     A = VA(IA)
     B = VB(IA)
     C = VC(IA)
     I = I_BAK
     J = J_BAK
     K = K_BAK
     X = X_BAK
     Y = Y_BAK
     Z = Z_BAK

     !+-----------------------------------------------------------------
     ! If the scattered photon travels towards the sensor, then:
     !+-----------------------------------------------------------------
     !IF (((k.LE.k_detector).AND.(vc(ia).GT.0.D0)).OR.((k.GE.k_detector).AND.(vc(ia).LT.0.D0))) THEN
        !IF (CHATTY) PRINT*,'Scattered photon travels towards sensor'

        IF (S_PROB(IA,1) .GT. EPS) THEN

           IF (CHATTY) PRINT*,'S_PROB(IA,1)=',S_PROB(IA,1)
           IF (CHATTY) PRINT*,'(S_PROB(IA,1) .GT. EPS)'

           !+-----------------------------------------------------------------
           ! If the sign of direction changed, simply reverse x, y and/or z
           !+-----------------------------------------------------------------
           IF (A*A_BAK .LT. 0.D0) X = MAX(0.D0, XUNIT - X)
           IF (B*B_BAK .LT. 0.D0) Y = MAX(0.D0, YUNIT - Y)
           IF (C*C_BAK .LT. 0.D0) Z = MAX(0.D0, ZUNIT(K) - Z)

           IF (CHATTY) PRINT*,'X=',X,'Y=',Y,'Z=',Z

           !+-----------------------------------------------------------------
           ! Zero TAU_V and path
           !+-----------------------------------------------------------------
           TAU_V      = 0.D0

           !+-----------------------------------------------------------------
           ! Get signs, abs values and inverses
           !+-----------------------------------------------------------------
           I_SIGN     = INT(SIGN(1.D0, A))
           J_SIGN     = INT(SIGN(1.D0, B))
           K_SIGN     = INT(SIGN(1.D0, C))
           A_ABS      = ABS(A)
           B_ABS      = ABS(B)
           C_ABS      = ABS(C)
           A_INV      = 1.D0 / MAX(A_ABS, EPS)
           B_INV      = 1.D0 / MAX(B_ABS, EPS)
           C_INV      = 1.D0 / MAX(C_ABS, EPS)

           !+-----------------------------------------------------------------
           ! Store photon location, with signs in direction towards sensor
           !+-----------------------------------------------------------------           
           x_old = x
           y_old = y
           z_old = z

           !+-----------------------------------------------------------------
           ! Compute TAU_V from event to sensor
           !+----------------------------------------------------------------- 
           IF (CHATTY) PRINT*,'K=',K,'NZ=',NZ,'K_detector=',K_detector

           !As long as the photon is still in the scene
           !DO WHILE ((K .GE. 1 .AND. K .LE. NZ).AND.(k .NE. k_detector))
           DO WHILE (K .GE. 1 .AND. K .LT. NZ) !k = nz is vacuum layer

              !and if the photon is above the clouds             
              !IF (k.GE.lowest_k_of_only_gas) THEN

              !and if we do not save the path (only then we can do this)
              !IF (.NOT. save_path) THEN

              IF (precompute_detector_paths) THEN

!                !and if the path to the sensor is not blocked by aerosol            
                 IF (blocked_sensor_path(I,J,K,IA) .EQ. 0) THEN
                 !IF (.false.) THEN

                    !and if the photon is at a z-boundary
                    IF (DABS(Z-ZUNIT(K)).LE.EPS) THEN

                       !then we can just use TAU_gas
                       TAU_V = TAU_V + TAU_gas(K,IA)

                       k = k_detector !this line is not really needed

                       GOTO 7009

                    ENDIF
                 ENDIF

              ENDIF

              !Otherwise we compute the path optical thickness per cell
              BETA_VALUE = BETA(I,J,K)

              i_old = i
              j_old = j
              k_old = k

              CALL ADVANCE(NX, NY, NZ, I, J, K, X, Y, Z, &
                   &              I_SIGN, J_SIGN, K_SIGN, A_ABS, B_ABS, C_ABS, &
                   &              A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT, DIST)

              TAU_V = TAU_V + DIST * BETA_VALUE

7009       ENDDO

           !-----------------------------------------------------------------------------------
           ! Now see where this event would be registered!
           !-----------------------------------------------------------------------------------
           IF (CHATTY) PRINT*,'See where this event would be registered.'

           !In backward monte carlo mode, we register at the injection location
           IF (backward_mode) THEN
              i = i_inject
              j = j_inject
              GOTO 7010
           ENDIF

           !Check if the detector is in the same direction as the register grid, as seen from the event
           IF (k_register .EQ. k_detector) THEN
              direction = 1
           ELSE
              direction = -1
           ENDIF

           !Double calculation of the mapping for debugging
           IF (advance_test) THEN

              I = I_BAK
              J = J_BAK
              K = K_BAK

              IF (direction .EQ. -1) THEN
                 X = XUNIT-X_OLD
                 Y = YUNIT-Y_OLD
                 Z = ZUNIT(K)-Z_OLD
              ELSE
                 X = X_OLD
                 Y = Y_OLD
                 Z = Z_OLD
              ENDIF

              !DO WHILE ((K .GE. 1 .AND. K .LE. NZ).AND.(k .NE. k_detector))
              DO WHILE (K .GE. 1 .AND. K .LT. NZ) !k = nz is vacuum layer

                 CALL ADVANCE(NX, NY, NZ, I, J, K, X, Y, Z, &
                & direction*I_SIGN, direction*J_SIGN, direction*K_SIGN, A_ABS, B_ABS, C_ABS, &
                &              A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT, DIST)

              ENDDO

           ENDIF

           !Map the signal onto the nx*ny output register grid
           IF (K_BAK .GT. 1) THEN
              HEIGHT = SUM(ZUNIT(1:K_BAK-1))
           ELSE
              HEIGHT = 0.D0
           ENDIF

           IF (K_SIGN .GT. 0) THEN
              HEIGHT = HEIGHT+ZUNIT(K_BAK)-Z_OLD
           ELSE
              HEIGHT = HEIGHT+Z_OLD
           ENDIF

           IF (k_register .EQ. NZ) THEN
              HEIGHT = SUM(ZUNIT(1:NZ-1)) - HEIGHT
           ENDIF

           DELTA_X = A_ABS * HEIGHT / C_ABS
           DELTA_Y = B_ABS * HEIGHT / C_ABS

           IF (direction .EQ. -1) THEN
              DELTA_I = INT( ( X_OLD + DELTA_X ) / XUNIT )
              DELTA_J = INT( ( Y_OLD + DELTA_Y ) / YUNIT )
           ELSE
              DELTA_I = INT( ( (XUNIT-X_OLD) + DELTA_X ) / XUNIT )
              DELTA_J = INT( ( (YUNIT-Y_OLD) + DELTA_Y ) / YUNIT )
           ENDIF

           I_NEW = MODULO(I_BAK + direction*I_SIGN*DELTA_I + NX - 1, NX) + 1
           J_NEW = MODULO(J_BAK + direction*J_SIGN*DELTA_J + NY - 1, NY) + 1

           !Check if the calculation was correction
           IF (advance_test) THEN

           IF ((I_NEW .NE. I) .OR. (J_NEW .NE. J)) THEN
              IF ((I_NEW .NE. I) .AND. (MOD(ABS(X_OLD),XUNIT) .LT. eps)) THEN
                 PRINT *,'Warning: (I_NEW .NE. I): I_NEW=',I_NEW,'I=',I
                 PRINT *,'MOD(ABS(X_OLD),XUNIT)=',MOD(ABS(X_OLD),XUNIT)
              ELSEIF ((J_NEW .NE. J) .AND. (MOD(ABS(Y_OLD),YUNIT) .LT. eps)) THEN
                 PRINT *,'Warning: (J_NEW .NE. J): J_NEW=',J_NEW,'J=',J
                 PRINT *,'MOD(ABS(Y_OLD),YUNIT)=',MOD(ABS(Y_OLD),YUNIT)
              ELSE
                 PRINT *,'ERROR!'
                 PRINT *,'NX=',NX
                 PRINT *,'-1*I_SIGN*DELTA_I=',-1*I_SIGN*DELTA_I
                 PRINT *,'I_BAK + -1*I_SIGN*DELTA_I + NX - 1=',I_BAK - I_SIGN*DELTA_I + NX - 1
                 PRINT *,'MOD(I_BAK + -1*I_SIGN*DELTA_I + NX - 1, NX) + 1',MOD(I_BAK - I_SIGN*DELTA_I + NX - 1, NX) + 1
                 PRINT *,'MOD(-1,20)=',MODULO(-1,20)
                 PRINT *,'.........'
                 PRINT *,'X=',X,'Y=',Y,'Z=',Z
                 PRINT *,'HEIGHT=',HEIGHT
                 PRINT *,'A_ABS=',A_ABS
                 PRINT *,'B_ABS=',B_ABS
                 PRINT *,'C_ABS=',C_ABS
                 PRINT *,'XUNIT=',XUNIT
                 PRINT *,'YUNIT=',YUNIT
                 PRINT *,'X_OLD=',X_OLD
                 PRINT *,'Y_OLD=',Y_OLD
                 PRINT *,'Z_OLD=',Z_OLD
                 PRINT *,'(XUNIT-X_OLD)=',(XUNIT-X_OLD)
                 PRINT *,'(YUNIT-Y_OLD)=',(YUNIT-Y_OLD)
                 PRINT *,'X_BAK=',X_BAK
                 PRINT *,'Y_BAK=',Y_BAK
                 PRINT *,'DELTA_X=',DELTA_X
                 PRINT *,'DELTA_Y=',DELTA_Y
                 PRINT *,'I_BAK=',I_BAK
                 PRINT *,'J_BAK=',J_BAK
                 PRINT *,'K_BAK=',K_BAK
                 PRINT *,'DELTA_I=',DELTA_I
                 PRINT *,'DELTA_J=',DELTA_J
                 PRINT *,'I_SIGN=',I_SIGN
                 PRINT *,'J_SIGN=',J_SIGN
                 PRINT *,'K_SIGN=',K_SIGN
                 PRINT *,'direction=',direction
                 PRINT *,'SUM(ZUNIT(1:NZ-1))=',SUM(ZUNIT(1:NZ-1))
                 PRINT *,'SUM(ZUNIT(1:K_BAK-1))=',SUM(ZUNIT(1:K_BAK-1))
                 PRINT *,'k_register=',k_register
                 PRINT *,'k_detector=',k_detector
                 PRINT *,'NZ=',NZ
                 PRINT *,'I_NEW=',I_NEW,', I=',I
                 PRINT *,'J_NEW=',J_NEW,', J=',J
                 READ(*,*)
              ENDIF
           ENDIF
           ENDIF

           !Update I and J
           I = I_NEW
           J = J_NEW

           !+-----------------------------------------------------------------
           ! Compute contribution
           !+-----------------------------------------------------------------
7010       IF (CHATTY) PRINT*,'Compute contribution'

           CONTRIB_STOKES(:) = S_PROB(IA,:) * EXP(-TAU_V)

           !+-----------------------------------------------------------------
           ! Store the origins of the single and multiple scattering contribtions
           !+-----------------------------------------------------------------           
           CONTRIB_STOKES_3D(:,i_store,j_store,k_store,ia) = &
                       & CONTRIB_STOKES_3D(:,i_store,j_store,k_store,ia) + CONTRIB_STOKES(:)
           IF (nscat .EQ. 1) THEN
              SINGSCAT_CONTRIB_STOKES_3D(:,i_store,j_store,k_store,ia) = &
                       & SINGSCAT_CONTRIB_STOKES_3D(:,i_store,j_store,k_store,ia) + CONTRIB_STOKES(:)
           ELSEIF (nscat .GT. 1) THEN
              MULTSCAT_CONTRIB_STOKES_3D(:,i_store,j_store,k_store,ia) = &
                       & MULTSCAT_CONTRIB_STOKES_3D(:,i_store,j_store,k_store,ia) + CONTRIB_STOKES(:)
           ELSE
              PRINT*,'ERROR: nscat should be larger than 0 when computing contribution'
              STOP
           ENDIF

           !+-----------------------------------------------------------------
           ! Compute radiances and Stokes
           !+-----------------------------------------------------------------              

           IF (chatty) PRINT*,'STOKES(:,I,J,IA)=',STOKES(:,I,J,IA),'CONTRIB_STOKES(:)',CONTRIB_STOKES(:)

           STOKES(:,I,J,IA) = STOKES(:,I,J,IA) + CONTRIB_STOKES(:)

           IF (nscat .EQ. 1) THEN
              STOKES_SING(:,I,J,IA)  = STOKES_SING(:,I,J,IA) + CONTRIB_STOKES(:)
           ELSE IF (nscat .GT. 1) THEN
              STOKES_MULT(:,I,J,IA)  = STOKES_MULT(:,I,J,IA) + CONTRIB_STOKES(:)
           ENDIF

           IF (chatty) PRINT*,'STOKES(:,I,J,IA)=',STOKES(:,I,J,IA),'CONTRIB_STOKES(:)',CONTRIB_STOKES(:)
           IF (pauses) READ(*,*)

        ENDIF

  ENDDO

  !+-----------------------------------------------------------------
  ! Restore actual photon position and direction.
  !+-----------------------------------------------------------------
  A = A_BAK
  B = B_BAK
  C = C_BAK
  I = I_BAK
  J = J_BAK
  K = K_BAK
  X = X_BAK
  Y = Y_BAK
  Z = Z_BAK

  RETURN

  END SUBROUTINE CALC_RADIANCES_AND_SCATS



!+-----------------------------------------------------------------
! SUBROUTINE calc_radiances_and_paths_pol
!+-----------------------------------------------------------------
! Compute reflectance contribution of this scattering event by
! by propagating the scattered light towards the detector.
! Version that saves air mass factors.
!
! Victor Trees
!+-----------------------------------------------------------------
SUBROUTINE calc_radiances_and_paths_pol(nx, ny, nz, nmat, n_angles, &
     &     nx_2, ny_2, nz_2, nmat_2, n_angles_2, &
     &     a, b, c, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
     &     TAU_gas,lowest_k_of_only_gas, &
     &     blocked_sensor_path, precompute_detector_paths,&
     &     va, vb, vc ,k_detector,&
     &     S_prob,Stokes,path_layers, &
     !&     nz_box,weight,path_boxes,     &
     &     weighted_path_store,&!weighted_boxpath_store,
     &     sum_of_weights_store,k_register,backward_mode,i_inject,j_inject)

  IMPLICIT NONE

  ! Input -------------------------------------------------------------------------------------------
  INTEGER,INTENT(IN) :: nx,ny,nz              !Number of grid cells in x, y and z directions (-)
  INTEGER,INTENT(IN) :: nmat                  !Size of the Stokes vector: 1, 3 or 4 (-)
  INTEGER,INTENT(IN) :: n_angles              !Number of viewing zenith angles (-)
  INTEGER,INTENT(IN) :: nx_2,ny_2,nz_2        !Similar, but only relevant if save_path=.true.
  INTEGER,INTENT(IN) :: nmat_2                !Similar, but only relevant if save_path=.true.
  INTEGER,INTENT(IN) :: n_angles_2            !Similar, but only relevant if save_path=.true.
  INTEGER,INTENT(IN) :: k_detector            !Layer index of the detector location (-)
  INTEGER,INTENT(IN) :: k_register            !Layer index of the save grid (-)
  INTEGER,INTENT(IN) :: lowest_k_of_only_gas  !Lowest layer index above which only gas occurs (-)
  INTEGER,INTENT(IN) :: i_inject,j_inject     !Column indices of photon injection in backward mode (-)
  !INTEGER,INTENT(IN) :: nz_box               !Number of layers for which box AMF are saved (-)      
  INTEGER,INTENT(IN) :: blocked_sensor_path(nx,ny,nz,n_angles) !Indicates if path to detector is blocked by aerosol (-)

  REAL*8,INTENT(IN) :: BETA(NX, NY, NZ)        !Extinction coefficient (1/km)
  REAL*8,INTENT(IN) :: xunit, yunit, ZUNIT(NZ) !Dimensions of the grid cells in the x, y and z diretion (km)
  REAL*8,INTENT(IN) :: VA(N_ANGLES), VB(N_ANGLES), VC(N_ANGLES) !x, y and z direction towards detector (-)
  REAL*8,INTENT(IN) :: S_PROB(N_ANGLES,nmat)   !Reflected Stokes vector at location of event (-)
  REAL*8,INTENT(IN) :: TAU_gas(NZ,N_ANGLES)    !Optical depth of path through gas towards detector (-)
  REAL*8,INTENT(IN) :: path_layers(NZ_2,N_ANGLES_2) !Accumulated paths inside the layers (km)
  !REAL*8,INTENT(IN) :: path_boxes(NX,NY,NZ_BOX,N_ANGLES) !Accumulated paths inside the cells (km)

  LOGICAL,INTENT(IN) :: backward_mode          !Indicates if the code runs in backward mode
  LOGICAL,INTENT(IN) :: precompute_detector_paths !Indicate if detector paths are precomputed

  ! In- and output ---------------------------------------------------------------------------------
  INTEGER,INTENT(INOUT) :: i,j,k               !Indices of the x, y and z grid cells

  REAL*8,INTENT(INOUT) :: a,b,c                !x, y and z directions of the photon before scattering (-)
  REAL*8,INTENT(INOUT) :: x,y,z                !Straight distances from photon to next x, y and z boundaries (km)
  REAL*8,INTENT(INOUT) :: STOKES(nmat,NX, NY, N_ANGLES) !Accumulated signal to the detector (-)
  REAL*8,INTENT(INOUT) :: weighted_path_store(NX_2,NY_2,NZ_2,N_ANGLES_2) !Accumulated weighted path through cells (km)
  !REAL*8,INTENT(INOUT) :: weighted_boxpath_store(NX,NY,NX,NY,NZ_BOX,N_ANGLES) !Accumulated weighted path through cells per cell (km)
  REAL*8,INTENT(INOUT) :: sum_of_weights_store(nx_2,ny_2,n_angles_2) !Accumulated signal through cells (-)

  ! Internal --------------------------------------------------------------------------------------
  INTEGER :: DELTA_I,DELTA_J,I_NEW,J_NEW,direction
  INTEGER :: i_bak,j_bak,k_bak,ia
  INTEGER :: i_sign,j_sign,k_sign
  INTEGER :: I_OLD,J_OLD,K_OLD

  REAL*8 :: tau_v,dist, beta_value
  REAL*8 :: a_bak,b_bak,c_bak
  REAL*8 :: a_abs,b_abs,c_abs,a_inv,b_inv,c_inv
  REAL*8 :: x_bak,y_bak,z_bak
  REAL*8 :: HEIGHT,DELTA_X,DELTA_Y,x_old,y_old,z_old
  REAL*8 :: CONTRIB_STOKES(nmat),CONTRIB_STOKES_NORM(nmat)
  REAL*8  :: path_layers_angle(nz_2)
  !REAL*8  :: path_boxes_angle(nx,ny,nz_box)

  REAL*8,PARAMETER :: eps = 1.0D-16

  LOGICAL :: chatty,pauses,advance_test

  chatty = .false.      !Print comments to screen
  pauses = .false.      !Pause script regularly
  advance_test = .false.!Double check if projection is okay

  !+-----------------------------------------------------------------
  ! Record actual photon position and direction.
  !+-----------------------------------------------------------------
  A_BAK = A
  B_BAK = B
  C_BAK = C
  I_BAK = I
  J_BAK = J
  K_BAK = K
  X_BAK = X
  Y_BAK = Y
  Z_BAK = Z

  !+-----------------------------------------------------------------
  ! Loop over each viewing direction.  Reset the photon position and 
  ! direction.  Zero TAU_V, the optical thickness remaining in the 
  ! model in the given direction.
  !+-----------------------------------------------------------------
  DO IA=1,N_ANGLES
     IF (CHATTY) PRINT*,'IA=',IA

     A = VA(IA)
     B = VB(IA)
     C = VC(IA)
     I = I_BAK
     J = J_BAK
     K = K_BAK
     X = X_BAK
     Y = Y_BAK
     Z = Z_BAK

     !+-----------------------------------------------------------------
     ! If the scattered photon travels towards the sensor, then:
     !+-----------------------------------------------------------------
     !IF (((k.LE.k_detector).AND.(vc(ia).GT.0.D0)).OR.((k.GE.k_detector).AND.(vc(ia).LT.0.D0))) THEN
        !IF (CHATTY) PRINT*,'Scattered photon travels towards sensor'

        IF (S_PROB(IA,1) .GT. EPS) THEN

           IF (CHATTY) PRINT*,'S_PROB(IA,1)=',S_PROB(IA,1)
           IF (CHATTY) PRINT*,'(S_PROB(IA,1) .GT. EPS)'

           !+-----------------------------------------------------------------
           ! If the sign of direction changed, simply reverse x, y and/or z
           !+-----------------------------------------------------------------
           IF (A*A_BAK .LT. 0.D0) X = MAX(0.D0, XUNIT - X)
           IF (B*B_BAK .LT. 0.D0) Y = MAX(0.D0, YUNIT - Y)
           IF (C*C_BAK .LT. 0.D0) Z = MAX(0.D0, ZUNIT(K) - Z)

           IF (CHATTY) PRINT*,'X=',X,'Y=',Y,'Z=',Z

           !+-----------------------------------------------------------------
           ! Zero TAU_V and path
           !+-----------------------------------------------------------------
           TAU_V      = 0.D0

           path_layers_angle = 0.D0
           !path_boxes_angle= 0.D0

           !+-----------------------------------------------------------------
           ! Get signs, abs values and inverses
           !+-----------------------------------------------------------------
           I_SIGN     = INT(SIGN(1.D0, A))
           J_SIGN     = INT(SIGN(1.D0, B))
           K_SIGN     = INT(SIGN(1.D0, C))
           A_ABS      = ABS(A)
           B_ABS      = ABS(B)
           C_ABS      = ABS(C)
           A_INV      = 1.D0 / MAX(A_ABS, EPS)
           B_INV      = 1.D0 / MAX(B_ABS, EPS)
           C_INV      = 1.D0 / MAX(C_ABS, EPS)

           !+-----------------------------------------------------------------
           ! Store photon location, with signs in direction towards sensor
           !+-----------------------------------------------------------------           
           x_old = x
           y_old = y
           z_old = z

           !+-----------------------------------------------------------------
           ! Compute TAU_V from event to sensor
           !+----------------------------------------------------------------- 
           IF (CHATTY) PRINT*,'K=',K,'NZ=',NZ,'K_detector=',K_detector

           !As long as the photon is still in the scene
           !DO WHILE ((K .GE. 1 .AND. K .LE. NZ).AND.(k .NE. k_detector))
           DO WHILE (K .GE. 1 .AND. K .LT. NZ)

              !We compute the path optical thickness per cell
              BETA_VALUE = BETA(I,J,K)

              i_old = i
              j_old = j
              k_old = k

              CALL ADVANCE(NX, NY, NZ, I, J, K, X, Y, Z, &
                   &              I_SIGN, J_SIGN, K_SIGN, A_ABS, B_ABS, C_ABS, &
                   &              A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT, DIST)

              TAU_V = TAU_V + DIST * BETA_VALUE

              ! Save the light path for this layer
              path_layers_angle(k_old) = path_layers_angle(k_old) + DIST
              !IF (k_old .LE. nz_box) THEN
              !   path_boxes_angle(i_old,j_old,k_old) = path_boxes_angle(i_old,j_old,k_old) + DIST
              !ENDIF

           ENDDO

           !-----------------------------------------------------------------------------------
           ! Now see where this event would be registered!
           !-----------------------------------------------------------------------------------
           IF (CHATTY) PRINT*,'See where this event would be registered.'

           !In backward monte carlo mode, we register at the injection location
           IF (backward_mode) THEN
              i = i_inject
              j = j_inject
              GOTO 6010
           ENDIF

           !Check if the detector is in the same direction as the register grid, as seen from the event
           IF (k_register .EQ. k_detector) THEN
              direction = 1
           ELSE
              direction = -1
           ENDIF

           !Double calculation of the mapping for debugging
           IF (advance_test) THEN

              I = I_BAK
              J = J_BAK
              K = K_BAK

              IF (direction .EQ. -1) THEN
                 X = XUNIT-X_OLD
                 Y = YUNIT-Y_OLD
                 Z = ZUNIT(K)-Z_OLD
              ELSE
                 X = X_OLD
                 Y = Y_OLD
                 Z = Z_OLD
              ENDIF

              !DO WHILE ((K .GE. 1 .AND. K .LE. NZ).AND.(k .NE. k_detector))
              DO WHILE (K .GE. 1 .AND. K .LT. NZ) !k = nz is vacuum layer

                 CALL ADVANCE(NX, NY, NZ, I, J, K, X, Y, Z, &
                & direction*I_SIGN, direction*J_SIGN, direction*K_SIGN, A_ABS, B_ABS, C_ABS, &
                &              A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT, DIST)

              ENDDO

           ENDIF

           !Map the signal onto the nx*ny output register grid
           IF (K_BAK .GT. 1) THEN
              HEIGHT = SUM(ZUNIT(1:K_BAK-1))
           ELSE
              HEIGHT = 0.D0
           ENDIF

           IF (K_SIGN .GT. 0) THEN
              HEIGHT = HEIGHT+ZUNIT(K_BAK)-Z_OLD
           ELSE
              HEIGHT = HEIGHT+Z_OLD
           ENDIF

           IF (k_register .EQ. NZ) THEN
              HEIGHT = SUM(ZUNIT(1:NZ-1)) - HEIGHT
           ENDIF

           DELTA_X = A_ABS * HEIGHT / C_ABS
           DELTA_Y = B_ABS * HEIGHT / C_ABS

           IF (direction .EQ. -1) THEN
              DELTA_I = INT( ( X_OLD + DELTA_X ) / XUNIT )
              DELTA_J = INT( ( Y_OLD + DELTA_Y ) / YUNIT )
           ELSE
              DELTA_I = INT( ( (XUNIT-X_OLD) + DELTA_X ) / XUNIT )
              DELTA_J = INT( ( (YUNIT-Y_OLD) + DELTA_Y ) / YUNIT )
           ENDIF

           I_NEW = MODULO(I_BAK + direction*I_SIGN*DELTA_I + NX - 1, NX) + 1
           J_NEW = MODULO(J_BAK + direction*J_SIGN*DELTA_J + NY - 1, NY) + 1

           !Check if the calculation was correction
           IF (advance_test) THEN

           IF ((I_NEW .NE. I) .OR. (J_NEW .NE. J)) THEN
              IF ((I_NEW .NE. I) .AND. (MOD(ABS(X_OLD),XUNIT) .LT. eps)) THEN
                 PRINT *,'Warning: (I_NEW .NE. I): I_NEW=',I_NEW,'I=',I
                 PRINT *,'MOD(ABS(X_OLD),XUNIT)=',MOD(ABS(X_OLD),XUNIT)
              ELSEIF ((J_NEW .NE. J) .AND. (MOD(ABS(Y_OLD),YUNIT) .LT. eps)) THEN
                 PRINT *,'Warning: (J_NEW .NE. J): J_NEW=',J_NEW,'J=',J
                 PRINT *,'MOD(ABS(Y_OLD),YUNIT)=',MOD(ABS(Y_OLD),YUNIT)
              ELSE
                 PRINT *,'ERROR!'
                 PRINT *,'NX=',NX
                 PRINT *,'-1*I_SIGN*DELTA_I=',-1*I_SIGN*DELTA_I
                 PRINT *,'I_BAK + -1*I_SIGN*DELTA_I + NX - 1=',I_BAK - I_SIGN*DELTA_I + NX - 1
                 PRINT *,'MOD(I_BAK + -1*I_SIGN*DELTA_I + NX - 1, NX) + 1',MOD(I_BAK - I_SIGN*DELTA_I + NX - 1, NX) + 1
                 PRINT *,'MOD(-1,20)=',MODULO(-1,20)
                 PRINT *,'.........'
                 PRINT *,'X=',X,'Y=',Y,'Z=',Z
                 PRINT *,'HEIGHT=',HEIGHT
                 PRINT *,'A_ABS=',A_ABS
                 PRINT *,'B_ABS=',B_ABS
                 PRINT *,'C_ABS=',C_ABS
                 PRINT *,'XUNIT=',XUNIT
                 PRINT *,'YUNIT=',YUNIT
                 PRINT *,'X_OLD=',X_OLD
                 PRINT *,'Y_OLD=',Y_OLD
                 PRINT *,'Z_OLD=',Z_OLD
                 PRINT *,'(XUNIT-X_OLD)=',(XUNIT-X_OLD)
                 PRINT *,'(YUNIT-Y_OLD)=',(YUNIT-Y_OLD)
                 PRINT *,'X_BAK=',X_BAK
                 PRINT *,'Y_BAK=',Y_BAK
                 PRINT *,'DELTA_X=',DELTA_X
                 PRINT *,'DELTA_Y=',DELTA_Y
                 PRINT *,'I_BAK=',I_BAK
                 PRINT *,'J_BAK=',J_BAK
                 PRINT *,'K_BAK=',K_BAK
                 PRINT *,'DELTA_I=',DELTA_I
                 PRINT *,'DELTA_J=',DELTA_J
                 PRINT *,'I_SIGN=',I_SIGN
                 PRINT *,'J_SIGN=',J_SIGN
                 PRINT *,'K_SIGN=',K_SIGN
                 PRINT *,'direction=',direction
                 PRINT *,'SUM(ZUNIT(1:NZ-1))=',SUM(ZUNIT(1:NZ-1))
                 PRINT *,'SUM(ZUNIT(1:K_BAK-1))=',SUM(ZUNIT(1:K_BAK-1))
                 PRINT *,'k_register=',k_register
                 PRINT *,'k_detector=',k_detector
                 PRINT *,'NZ=',NZ
                 PRINT *,'I_NEW=',I_NEW,', I=',I
                 PRINT *,'J_NEW=',J_NEW,', J=',J
                 READ(*,*)
              ENDIF
           ENDIF
           ENDIF

           !Update I and J
           I = I_NEW
           J = J_NEW

           !+-----------------------------------------------------------------
           ! Compute contribution
           !+-----------------------------------------------------------------
6010       IF (CHATTY) PRINT*,'Compute contribution'

           !IF (wrap==0) THEN

           CONTRIB_STOKES(:) = S_PROB(IA,:) * EXP(-TAU_V)

           !+-----------------------------------------------------------------
           ! Store the light paths per layer (k) for this pixel (i,j,ia)
           !+-----------------------------------------------------------------
           weighted_path_store(i,j,:,ia) = weighted_path_store(i,j,:,ia) + &
                                          & CONTRIB_STOKES(1) * (path_layers(:,ia) + path_layers_angle(:))

           sum_of_weights_store(i,j,ia) = sum_of_weights_store(i,j,ia) + CONTRIB_STOKES(1)


           !This line is very time consuming:
           !weighted_boxpath_store(i,j,:,:,:,ia) = weighted_boxpath_store(i,j,:,:,:,ia) + &
           !                                  & CONTRIB_STOKES(1) * (path_boxes(:,:,:,ia) + path_boxes_angle(:,:,:))

           !+-----------------------------------------------------------------
           ! Compute radiances and Stokes
           !+-----------------------------------------------------------------              

           IF (chatty) PRINT*,'STOKES(:,I,J,IA)=',STOKES(:,I,J,IA),'CONTRIB_STOKES(:)',CONTRIB_STOKES(:)

           STOKES(:,I,J,IA) = STOKES(:,I,J,IA) + CONTRIB_STOKES(:)

           !IF (save_path) THEN
           !   weighted_path_store(i,j,:,ia) = weighted_path_store(i,j,:,ia) + CONTRIB_STOKES(1)*path_layers(:,ia)
           !ENDIF

           IF (chatty) PRINT*,'STOKES(:,I,J,IA)=',STOKES(:,I,J,IA),'CONTRIB_STOKES(:)',CONTRIB_STOKES(:)
           IF (pauses) READ(*,*)

        ENDIF

  ENDDO

  !+-----------------------------------------------------------------
  ! Restore actual photon position and direction.
  !+-----------------------------------------------------------------
  A = A_BAK
  B = B_BAK
  C = C_BAK
  I = I_BAK
  J = J_BAK
  K = K_BAK
  X = X_BAK
  Y = Y_BAK
  Z = Z_BAK

  RETURN

  END SUBROUTINE CALC_RADIANCES_AND_PATHS_POL


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
     &     A_INV, B_INV, C_INV, XUNIT, YUNIT, ZUNIT, DIST)

  IMPLICIT NONE

  ! Input --------------------------------------------------------
  INTEGER,INTENT(IN)     :: NX,NY,NZ,I_SIGN,J_SIGN,K_SIGN
  REAL*8,INTENT(IN)      :: A_ABS,B_ABS,C_ABS
  REAL*8,INTENT(IN)      :: A_INV,B_INV,C_INV,XUNIT,YUNIT,ZUNIT(NZ)

  ! In- and output -----------------------------------------------
  INTEGER,INTENT(INOUT)  :: I,J,K
  REAL*8,INTENT(INOUT)   :: X,Y,Z
  REAL*8,INTENT(INOUT)   :: DIST
 
  ! Internal -----------------------------------------------------
  REAL*8                 :: DIST_X,DIST_Y,DIST_Z

  !---------------------------------------------------------------
  ! Get distances to x, y and z boundaries along photon path
  ! and compute the minimum
  !---------------------------------------------------------------
  DIST_X = X * A_INV
  DIST_Y = Y * B_INV
  DIST_Z = Z * C_INV
  DIST   = MIN(DIST_X, DIST_Y, DIST_Z)

  !---------------------------------------------------------------
  ! Determine intersected boundary and update photon location
  !---------------------------------------------------------------
  IF (DIST_X .EQ. DIST) THEN ! the photon reaches the x boundary
     I = MOD(I + I_SIGN + NX - 1, NX) + 1
     X = XUNIT
     Y = MAX(0.D0, Y - B_ABS * DIST)
     Z = MAX(0.D0, Z - C_ABS * DIST)
  ELSE IF (DIST_Y .EQ. DIST) THEN ! the photon reaches the y boundary
     J = MOD(J + J_SIGN + NY - 1, NY) + 1
     X = MAX(0.D0, X - A_ABS * DIST)
     Y = YUNIT
     Z = MAX(0.D0, Z - C_ABS * DIST)
  ELSE                      ! the photon reaches the z boundary
     K = K + K_SIGN
     X = MAX(0.D0, X - A_ABS * DIST)
     Y = MAX(0.D0, Y - B_ABS * DIST)
     IF (K .GE. 1 .AND. K .LE. NZ) Z = ZUNIT(K)
  ENDIF

  RETURN
END SUBROUTINE ADVANCE


END MODULE mc_routines_cartesian



