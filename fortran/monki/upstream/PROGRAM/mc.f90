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
! SUBROUTINE monte_carlo
!
! Monte Carlo subroutine of MONKI.
!
! Author and lead developer: Victor J.H. Trees
!-----------------------------------------------------------------------
SUBROUTINE monte_carlo(nx, ny, nz, nmat, n_angles,& 
                & nx_1, ny_1, nz_1, nmat_1, n_angles_1,&
                & nx_2, ny_2, nz_2, nmat_2, n_angles_2,&
                & n_types, new_scatt_data_pol, &
                & n_photons, &
                & vartheta0,phi0, &
                & vartheta,phi, &
                & xunit, yunit, zunit, &
                & surface,surface_alb,i_surface_elev, &
                & beta, ext_by_species, ssalb_gas, depol, &
                & k_detector, &
                & TAU_gas,lowest_k_of_only_gas, &
                & blocked_sensor_path,precompute_detector_paths, &
                & Stokes, &
                & save_path,&!nz_box,&
                & weighted_path_store,&
                !& weighted_boxpath_store,&
                & sum_of_weights_store,&
                & save_scats,&
                & CONTRIB_STOKES_3D,&
                & SINGSCAT_CONTRIB_STOKES_3D,&
                & MULTSCAT_CONTRIB_STOKES_3D,&
                & Stokes_sing,Stokes_mult,&
                & k_register,k_source,backward_mode)

  !---------------------------------------------------------------------
  ! Import modules
  !---------------------------------------------------------------------        
  USE mc_routines_general
  USE mc_routines_cartesian
  USE new_data_types
  USE physical_parameters
  
  !---------------------------------------------------------------------
  ! Declare variables
  !---------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Input --------------------------------------------------------------
  INTEGER,INTENT(IN) :: nx,ny,nz              !Number of grid cells in x, y and z direction (-)
  INTEGER,INTENT(IN) :: nmat                  !Size of the Stokes vector: 1, 3 or 4 (-)
  INTEGER,INTENT(IN) :: n_angles              !Number of viewing zentih angles (-)
  INTEGER,INTENT(IN) :: nx_1,ny_1,nz_1        !Similar, but only relevant if save_scats = .true.
  INTEGER,INTENT(IN) :: nmat_1                !Similar, but only relevant if save_scats = .true.
  INTEGER,INTENT(IN) :: n_angles_1            !Similar, but only relevant if save_scats = .true.
  INTEGER,INTENT(IN) :: nx_2,ny_2,nz_2        !Similar, but only relevant if save_path = .true.
  INTEGER,INTENT(IN) :: nmat_2                !Similar, but only relevant if save_path = .true.
  INTEGER,INTENT(IN) :: n_angles_2            !Similar, but only relevant if save_path = .true.
  INTEGER,INTENT(IN) :: n_types               !Number of cloud droplet types (-)
  INTEGER,INTENT(IN) :: k_register            !Layer index where photons are registered (-)
  INTEGER,INTENT(IN) :: k_source              !Layer index of photon source (-)
  INTEGER,INTENT(IN) :: k_detector            !Layer index of detector (-)
  INTEGER,INTENT(IN) :: lowest_k_of_only_gas  !Layer index above which only gas occurs (-)
  INTEGER,INTENT(IN) :: blocked_sensor_path(nx,ny,nz,n_angles) !Flags (1) grid cells whose paths towards detector are blocked by aerosols (-)
  INTEGER,INTENT(IN) :: surface(nx, ny)       !Surface type indicator, -1 is Lambertian (-)
  INTEGER,INTENT(IN) :: i_surface_elev(nx, ny)!Surface elevation matrix (-)
  !INTEGER,INTENT(IN) :: nz_box               !Number of layers for which box AMF are saved

  INTEGER*8,INTENT(IN) :: n_photons             !Number of photons for this round and thread (-)

  REAL*8,INTENT(IN) :: xunit,yunit,zunit(nz)  !Dimensions of the grid cells in x, y and z direction (km)
  REAL*8,INTENT(IN) :: vartheta0              !Solar zenith angle (deg)
  REAL*8,INTENT(IN) :: vartheta(n_angles)     !Viewing zenith angles (deg)
  REAL*8,INTENT(IN) :: phi0                   !Solar azimuth angle (deg)
  REAL*8,INTENT(IN) :: phi                    !Viewing azimuth angle (deg)
  REAL*8,INTENT(IN) :: TAU_gas(nz,n_angles)   !Optical depth of the gas paths towards detector (-)
  REAL*8,INTENT(IN) :: surface_alb(nx, ny)    !Surface albedo (-)
  REAL*8,INTENT(IN) :: beta(nx, ny, nz)       !Extinction coefficient (1/km)
  REAL*8,INTENT(IN) :: ssalb_gas(nz)          !Single scattering albedo of the gas (-)
  REAL*8,INTENT(IN) :: depol                  !Depolarization factor (-)
  REAL*8,INTENT(IN) :: ext_by_species(n_types+1, nx, ny, nz) !Extinction coefficient per species, normalized by total Beta (-)

  LOGICAL,INTENT(IN) :: precompute_detector_paths !Logical to indicate whether to precompute detector paths
  LOGICAL,INTENT(IN) :: save_path             !Logical to indicate whether to save the AMFs
  LOGICAL,INTENT(IN) :: save_scats            !Logical to indicate whether to save single and multiple scat seperately
  LOGICAL,INTENT(IN) :: backward_mode         !Logical to indicate whether code runs in backward mode

  TYPE(new_scatterer_info_pol),INTENT(IN)  :: new_scatt_data_pol(n_types) !Data type containing droplet single scattering properties

  ! Output -------------------------------------------------------------
  REAL*8,INTENT(OUT) :: Stokes(nmat,nx,ny,n_angles)  !Accumulated reflectance or transmittance vector (-)
  REAL*8,INTENT(OUT) :: CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1) !3D contribution of signal (-)
  REAL*8,INTENT(OUT) :: SINGSCAT_CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1) !3D single scattering contribution to signal (-)
  REAL*8,INTENT(OUT) :: MULTSCAT_CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1) !3D multiple scattering contribution to signal
  REAL*8,INTENT(OUT) :: Stokes_sing(nmat_1,nx_1,ny_1,n_angles_1) !Single scattering contribution to signal
  REAL*8,INTENT(OUT) :: Stokes_mult(nmat_1,nx_1,ny_1,n_angles_1) !Multiple scattering contribution to signal
  REAL*8,INTENT(OUT) :: weighted_path_store(nx_2,ny_2,nz_2,n_angles_2) !Accumulated light path through cells (km)
  !REAL*8,INTENT(OUT) :: weighted_boxpath_store(nx_2,ny_2,nx_2,ny_2,nz_box_2,n_angles_2) !Accumulated light path through cells per cell (km)
  REAL*8,INTENT(OUT) :: sum_of_weights_store(nx_2,ny_2,n_angles_2) !Accumulated signal through cells (-) 

  ! Internal -----------------------------------------------------------
  INTEGER :: i_photon,i_column,nscat,i_surf_elev
  INTEGER :: ia,i,j,k,i_sign,j_sign,k_sign,is,i_old,j_old,k_old
  INTEGER :: i_inject,j_inject,i_store_singscat,j_store_singscat,k_store_singscat

  INTEGER*8 :: n_columns,n_distribute

  REAL :: work,rtemp

  REAL*8 :: S(nmat),S_back(nmat),S_prob(n_angles,nmat),x_comp,y_comp,z_comp
  REAL*8 :: total_S1,total_S1_sing,total_S1_mult
  REAL*8 :: P(nmat,nmat),P_surf(nmat,nmat),L1(nmat,nmat),L2(nmat,nmat),F_sca(nmat,nmat)
  REAL*8 :: cos_theta_sca,dphi_rad,phi_view_rad,phi_in_rad,phi_out_rad
  REAL*8 :: cos_theta_in,cos_theta_out,phi0_rad,phi_rad,vartheta0_rad
  REAL*8 :: vartheta_rad(n_angles)
  REAL*8 :: a_bak,b_bak,c_bak,a0,b0,c0,a,b,c,a_abs,b_abs,c_abs
  REAL*8 :: a_inv,b_inv,c_inv,a_int,b_int,c_int,weight,remains
  REAL*8 :: va(n_angles),vb(n_angles),vc(n_angles)
  REAL*8 :: beta_value,dist_event,dist
  REAL*8 :: wp,x,y,z,mu0mu_factor
  REAL*8 :: path_layers(nz_2,n_angles_2)
  REAL*8 :: S0(nmat),weight_matrix(nmat,nmat),As(nmat,nmat)
  !REAL*8 :: path_boxes(nx,ny,nz_box,n_angles)
  REAL*8 :: stokes_criterium_error,zunit_surface

  REAL*8,PARAMETER  :: eps = 1.0D-16

  LOGICAL :: PDFerror,chatty,print_progress,pauses,checkinputs
  LOGICAL :: rough_surface

  !---------------------------------------------------------------------
  !************************ START CALCULATIONS *************************
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'in mc.f90, start calculations'

  !---------------------------------------------------------------------
  ! Debugging settings
  !---------------------------------------------------------------------
  pauses = .false.
  chatty = .false.
  print_progress = .false.
  checkinputs = .false.

  !---------------------------------------------------------------------
  ! Set allowed error in the Stokes criterium SQRT(Q**2+U**2+V**2).LE.I
  !---------------------------------------------------------------------  
  stokes_criterium_error = 1.D-2 !1.D-6 !1.D-8

  !---------------------------------------------------------------------
  ! Check inputs if needed
  !---------------------------------------------------------------------
  IF (checkinputs) THEN 
     PRINT *,'nx = ',nx
     PRINT *,'ny = ',ny
     PRINT *,'nz = ',nz
     PRINT *,'n_angles = ',n_angles
     PRINT *,'n_types = ',n_types
     READ(*,*)
     !READ(*,*)
     PRINT *,'n_photons = ',n_photons
     READ(*,*)
     PRINT *,'vartheta0 = ',vartheta0
     READ(*,*)
     PRINT *,'phi0 = ',phi0
     READ(*,*)
     PRINT *,'vartheta = ',vartheta
     READ(*,*)
     PRINT *,'phi = ',phi
     READ(*,*)
     PRINT *,'xunit = ',xunit
     PRINT *,'yunit = ',yunit
     PRINT *,'zunit = ',zunit
     READ(*,*)
     PRINT *,'surface = ',surface
     READ(*,*)
     PRINT *,'surface_alb = ',surface_alb
     READ(*,*)
     PRINT *,'i_surface_elev = ',i_surface_elev
     READ(*,*)
     PRINT *,'beta = ',beta
     READ(*,*)
     PRINT *,'ext_by_species = ',ext_by_species
     READ(*,*)
     PRINT *,'k_detector = ',k_detector
     READ(*,*)
     PRINT *,'Stokes = ',Stokes
     READ(*,*)
  ENDIF

  !---------------------------------------------------------------------
  ! Initialize variables
  !---------------------------------------------------------------------
  Stokes = 0.D0
  Stokes_sing = 0.D0
  Stokes_mult = 0.D0
  IF (save_path) THEN
     weighted_path_store = 0.D0
     !weighted_boxpath_store = 0.D0
     sum_of_weights_store = 0.D0
  ENDIF
  IF (save_scats) THEN
     contrib_stokes_3D = 0.D0     
     singscat_contrib_stokes_3D = 0.D0
     multscat_contrib_stokes_3D = 0.D0
  ENDIF

  !---------------------------------------------------------------------
  ! Compute number of columns
  !---------------------------------------------------------------------
  n_columns = nx * ny

  IF (chatty) PRINT *,'n_columns=',n_columns

  !---------------------------------------------------------------------
  ! Compute number of photons that can evenly be distributed among 
  ! the columns
  !---------------------------------------------------------------------
  n_distribute = INT(n_photons / n_columns) * n_columns

  IF (chatty) PRINT *,'n_distribute=',n_distribute  
 
  !---------------------------------------------------------------------
  ! Check if the surface is rough
  !---------------------------------------------------------------------
  IF (MAXVAL(i_surface_elev) .GT. 1) THEN
     rough_surface = .true.
  ELSE
     rough_surface = .false.
  ENDIF

  !---------------------------------------------------------------------
  ! Define a sursurface layer height were the light source is in 
  ! case it is located below the surface. Set it equal to the vacuum 
  ! layer height.
  !---------------------------------------------------------------------
  IF (k_source .EQ. 0) THEN
     zunit_surface = zunit(nz)
  ENDIF

  !---------------------------------------------------------------------
  ! Define the Lambertian surface reflectance matrix
  ! (The 1,1 element is just 1.D0, because the surface
  ! albedo is already adjusting the photon weight during
  ! surface interaction.)
  !---------------------------------------------------------------------  
  As = 0.D0
  As(1,1) = 1.D0

  !---------------------------------------------------------------------
  ! Convert angles from degrees to radians
  !---------------------------------------------------------------------   
  vartheta0_rad = DegtoRad*vartheta0
  phi0_rad = DegtoRad*phi0
  vartheta_rad(:) = DegtoRad*vartheta(:)
  phi_rad = DegtoRad*phi

  !---------------------------------------------------------------------
  ! Compute directional cosines of the solar angle
  !---------------------------------------------------------------------  
  a0 = DSIN(vartheta0_rad) * DCOS(phi0_rad)
  b0 = DSIN(vartheta0_rad) * DSIN(phi0_rad)
  c0 = DCOS(vartheta0_rad)

  !---------------------------------------------------------------------
  ! Set incident solar flux vector
  !---------------------------------------------------------------------   
  S0 = 0.D0
  S0(1) = 1.D0

  !---------------------------------------------------------------------
  ! Compute directional cosines of the viewing angles
  !---------------------------------------------------------------------   
  va(:)=DCOS(phi_rad)*DSIN(vartheta_rad(:))
  vb(:)=DSIN(phi_rad)*DSIN(vartheta_rad(:))
  vc(:)=DCOS(vartheta_rad(:))
  
  phi_view_rad = phi_rad

  !*********************************************************************
  ! +-------------------------------------------------------------------
  ! Start of simulation
  ! +-------------------------------------------------------------------
  ! Inject n_photons photons into the model one at a time.  
  ! Zero the scatter count and set the initial direction cosines.  
  ! Distribute the first n_distribute photons evenly among the top 
  ! cells, and randomize the starting cells of the remaining photons.  
  ! Randomize the starting positions within the cells.
  !
  ! (i,j,k) are the indices of the photon's current cell.  
  ! k = 1 at the bottom of the model.
  !
  ! (x,y,z) are the distances within the current cell to the
  ! boundaries being approached by the photon. For example, suppose
  ! the photon is 10 m from the top of a 50 m cell. If the photon is
  ! heading upward, z = 10 m. If it is heading downward, z = 40 m.
  !
  ! (a,b,c) are cosines that indicate the direction of the photon.
  !
  ! +-------------------------------------------------------------------
  !*********************************************************************

  !---------------------------------------------------------------------
  ! Start loop over photons
  !---------------------------------------------------------------------
  DO 4000 i_photon=1,n_photons   

     IF (chatty) PRINT *,''
     IF (chatty) PRINT *,'************************************'
     IF (chatty) PRINT *,'New photon injected at TOA'
     IF (chatty) PRINT *,'************************************'
  
     IF (chatty) PRINT *,'i_photon =',i_photon,' n_photons =',n_photons

     IF (print_progress) THEN
        IF (i_photon .EQ. INT(n_photons/10.0*0.01)) WRITE(*,'(A)',advance='no') ' 0.1 %'
        IF (i_photon .EQ. INT(n_photons/10.0*0.1)) WRITE(*,'(A)',advance='no') ' 1 %'     
        IF (i_photon .EQ. INT(n_photons/10.0*0.5)) WRITE(*,'(A)',advance='no') ' 5 %' 
        IF (i_photon .EQ. INT(n_photons/10.0*1.0)) WRITE(*,'(A)',advance='no') ' 10 %'
        IF (i_photon .EQ. INT(n_photons/10.0*1.5)) WRITE(*,'(A)',advance='no') ' 15 %'
        IF (i_photon .EQ. INT(n_photons/10.0*2.0)) WRITE(*,'(A)',advance='no') ' 20 %'
        IF (i_photon .EQ. INT(n_photons/10.0*2.5)) WRITE(*,'(A)',advance='no') ' 25 %'
        IF  (i_photon .EQ. INT(n_photons/10.0*3.0)) WRITE(*,'(A)',advance='no') ' 30 %'
        IF (i_photon .EQ. INT(n_photons/10.0*3.5)) WRITE(*,'(A)',advance='no') ' 35 %'
        IF (i_photon .EQ. INT(n_photons/10.0*4.0)) WRITE(*,'(A)',advance='no') ' 40 %'
        IF (i_photon .EQ. INT(n_photons/10.0*4.5)) WRITE(*,'(A)',advance='no') ' 45 %'
        IF (i_photon .EQ. INT(n_photons/10.0*5.0)) WRITE(*,'(A)',advance='no') ' 50 %'
        IF (i_photon .EQ. INT(n_photons/10.0*5.5)) WRITE(*,'(A)',advance='no') ' 55 %'
        IF (i_photon .EQ. INT(n_photons/10.0*6.0)) WRITE(*,'(A)',advance='no') ' 60 %'
        IF (i_photon .EQ. INT(n_photons/10.0*6.5)) WRITE(*,'(A)',advance='no') ' 65 %'
        IF (i_photon .EQ. INT(n_photons/10.0*7.0)) WRITE(*,'(A)',advance='no') ' 70 %'
        IF (i_photon .EQ. INT(n_photons/10.0*7.5)) WRITE(*,'(A)',advance='no') ' 75 %'
        IF (i_photon .EQ. INT(n_photons/10.0*8.0)) WRITE(*,'(A)',advance='no') ' 80 %'
        IF (i_photon .EQ. INT(n_photons/10.0*8.5)) WRITE(*,'(A)',advance='no') ' 85 %'
        IF (i_photon .EQ. INT(n_photons/10.0*9.0)) WRITE(*,'(A)',advance='no') ' 90 %'
        IF (i_photon .EQ. INT(n_photons/10.0*9.5)) WRITE(*,'(A)',advance='no') ' 95 %'
        IF (i_photon .EQ. INT(n_photons/10.0*10.0)) WRITE(*,'(A)',advance='no') ' 100 %'
     ENDIF     

     !---------------------------------------------------------------------
     ! Determine the index of the column the photon is injected into
     !---------------------------------------------------------------------
     IF (i_photon .LE. n_distribute) THEN

       !Distribute photons evenly among columns:      
       i_column = MOD(i_photon - 1, n_columns)  

     ELSE 

       !Distribute remaining photons randomly: 
       CALL RANDOM_NUMBER(work)
       i_column = INT(n_columns * work)   

     ENDIF    

     IF (chatty) PRINT *,'i_column=',i_column

     !---------------------------------------------------------------------
     ! Get 3D indices of the photon's current cell
     !--------------------------------------------------------------------- 
     i = MOD(i_column, nx)  + 1
     j = MIN(INT(i_column / nx) + 1,ny)
     
     !k = nz !always inject photon in vacuum layer above TOA
     k = k_source
     
     IF (chatty) PRINT *,'i=',i,'j=',j,'k=',k

     !save injection location in backward mode
     IF (backward_mode) THEN
        i_inject = i
        j_inject = j
     ENDIF

     !---------------------------------------------------------------------
     ! Set x, y and z:
     ! x and y are randomized
     ! z = zunit(k), unless k_sourse .EQ. 0 then photon starts below surface
     !---------------------------------------------------------------------
     CALL RANDOM_NUMBER(work)
     x = xunit * work
     CALL RANDOM_NUMBER(work)
     y = yunit * work

     IF (k_source .GT. 0) THEN 
        z = zunit(k)
     ELSE
        z = zunit_surface
     ENDIF

     IF (chatty) PRINT *,'x=',x,'y=',y,'z=',z
     
     !---------------------------------------------------------------------
     ! Initialize the Stokes vector for this photon
     !---------------------------------------------------------------------
     S = S0

     !---------------------------------------------------------------------
     ! Initialize the weight for this photon
     !---------------------------------------------------------------------
     weight = 1.D0

     weight_matrix = 0.D0
     weight_matrix(1,1)=1.D0
     IF (nmat .GT. 1) THEN
        weight_matrix(2,2)=1.D0
        weight_matrix(3,3)=1.D0
        IF (nmat .GT. 3) THEN
           weight_matrix(4,4)=1.D0
        ENDIF
     ENDIF

     !---------------------------------------------------------------------
     ! Initialize the light path for this photon
     !---------------------------------------------------------------------     
     IF (save_path) THEN
        path_layers = 0.D0
        !path_boxes = 0.D0
     ENDIF
     
     !---------------------------------------------------------------------
     ! Initialize the scattering counter and 1st order scattering location
     !--------------------------------------------------------------------- 
     IF (save_scats) THEN
        nscat = 0
        IF (backward_mode) THEN
           i_store_singscat = 0
           j_store_singscat = 0
           k_store_singscat = 0
        ENDIF
     ENDIF

     !---------------------------------------------------------------------
     ! Get solar angles vartheta0 and phi0 for the photon's column
     !---------------------------------------------------------------------
     !vartheta0 = DegtoRad*vartheta0_2d(i,j)
     !phi0 = DegtoRad*phi0_2d(i,j) 

     !IF (chatty) PRINT *,'vartheta0=',vartheta0,'phi0=',phi0

     !---------------------------------------------------------------------
     ! Compute directional cosines of the solar angles
     !---------------------------------------------------------------------  
     !a0 = DSIN(vartheta0) * DCOS(phi0)
     !b0 = DSIN(vartheta0) * DSIN(phi0)
     !c0 = DCOS(vartheta0)

     !---------------------------------------------------------------------
     ! Set the photon's direction equal to the direction of the sunlight
     !---------------------------------------------------------------------      
     a = a0
     b = b0
     c = c0
     phi_in_rad = phi0_rad

     IF (chatty) PRINT *,'a=',a,'b=',b,'c=',c

     ! ********************************************************************
     ! +-------------------------------------------------------------------
     ! photon advancing after injection at TOA or after an event
     ! +-------------------------------------------------------------------
     ! remains is -ln(RandomNumber)
     ! randomize its initial value.  
     ! return to 1000 after scattering events.
     ! +-------------------------------------------------------------------
     ! ********************************************************************

1000 CALL RANDOM_NUMBER(work)
     remains = -LOG(MAX(work,eps))

     IF (chatty) PRINT *,'remains=',remains

     ! ********************************************************************
     ! +-------------------------------------------------------------------
     ! Determine whether the photon will reach an event before a cell
     ! boundary. 
     ! If so, move the photons to where the event will occur,
     ! and goto 3000.  
     ! Return to 1010 if the photon changes direction, as after a 
     ! surface reflection.  
     ! Return to 1020 after the photon crosses a cell boundary (photon
     ! direction stays the same).
     ! +-------------------------------------------------------------------
     ! ********************************************************************
     
     !---------------------------------------------------------------------
     ! Compute signs, absolutes and inverses of the photon direction 
     !---------------------------------------------------------------------
1010 i_sign = INT(SIGN(1.D0, a))     !sign of a
     j_sign = INT(SIGN(1.D0, b))     !sign of b
     k_sign = INT(SIGN(1.D0, c))     !sign of c
     a_abs  = DABS(a)                !|a|
     b_abs  = DABS(b)                !|b|
     c_abs  = DABS(c)                !|c|
     a_inv  = 1.D0 / MAX(a_abs, eps) !1/|a|
     b_inv  = 1.D0 / MAX(b_abs, eps) !1/|b|
     c_inv  = 1.D0 / MAX(c_abs, eps) !1/|c|

1020 IF (chatty) PRINT *,''
     IF (chatty) PRINT *,'----------------------------------------'
     IF (chatty) PRINT *,'Check for atmospheric event in this cell'
     IF (chatty) PRINT *,'----------------------------------------'

     !--------------------------------------------------------------------
     ! Get extinction coefficient beta for this cell if we are not in 
     ! vacuum layer (nz) or below surface
     !--------------------------------------------------------------------     
     IF ((k.LT.nz) .AND. (k.GT.0)) THEN
        beta_value = beta(i,j,k)
     ELSE
        beta_value = 0.D0
     ENDIF

     IF (chatty) PRINT *,'beta_value=',beta_value

     IF (beta_value .GT. eps) THEN
       
        !-------------------------------------------------------------------- 
        ! Get photon path distance untill event
        !--------------------------------------------------------------------       
        dist_event = remains/beta_value

        IF (chatty) PRINT *,'dist_event=',dist_event 

        !--------------------------------------------------------------------
        ! IF dist_event*|a| < x 
        ! OR dist_event*|b| < y 
        ! OR dist_event*|c| < z
        ! THEN there is an atmospheric event
        !--------------------------------------------------------------------
        !IF (dist_event .LT. MIN(x*a_inv,y*b_inv,z*c_inv)) THEN  

        x_comp = x * a_inv
        y_comp = y * b_inv
        z_comp = z * c_inv
        
        IF (dist_event .LT. x_comp .AND. dist_event .LT. y_comp .AND. dist_event .LT. z_comp) THEN


           IF (chatty) PRINT *,'-----------------'      
           IF (chatty) PRINT *,'atmospheric event'
           IF (chatty) PRINT *,'-----------------' 

           !--------------------------------------------------------------------
           ! Move photon to position in cell where the event occurs
           !--------------------------------------------------------------------
           x = MAX(0.D0, x - a_abs * dist_event)
           y = MAX(0.D0, y - b_abs * dist_event)
           z = MAX(0.D0, z - c_abs * dist_event)
           
           IF (chatty) PRINT *,'x=',x,'y=',y,'z=',z

           !--------------------------------------------------------------------
           ! Save the light path for this layer
           !--------------------------------------------------------------------
           IF (save_path) THEN
                 path_layers(k,:) = path_layers(k,:) + dist_event
              !IF (k .LE. nz_box) THEN
              !   path_boxes(i,j,k,:) = path_boxes(i,j,k,:) + dist_event
              !ENDIF
           ENDIF           

           !--------------------------------------------------------------------
           ! Get cosines of viewing direction towards the detector for this
           ! new location of the photon
           !--------------------------------------------------------------------
           !va=DCOS(DegtoRad*phi)*DSIN(DegtoRad*vartheta)
           !vb=DSIN(DegtoRad*phi)*DSIN(DegtoRad*vartheta)
           !vc=DCOS(DegtoRad*vartheta)
  
           !IF (chatty) PRINT *,'va=',va,'vb=',vb,'vc=',vc

           !--------------------------------------------------------------------
           ! GOTO 3000: absorption or scattering event
           !--------------------------------------------------------------------
           GOTO 3000
        
        ELSE 
           IF (chatty) PRINT *,''     
           IF (chatty) PRINT *,'NO atmospheric event'
           IF (chatty) PRINT *,'' 
        ENDIF
        
     ELSE 
        IF (chatty) PRINT *,'beta_value < eps'  
     ENDIF
     
     IF (chatty) PRINT *,'----------------------------'   
     IF (chatty) PRINT *,'Photon reached cell boundary'
     IF (chatty) PRINT *,'----------------------------' 
     
     ! ******************************************************************
     ! +-----------------------------------------------------------------
     ! The photon will reach a cell boundary. 
     ! Advance the photon to the next cell and calculate dist, 
     ! the distance it travelled to get there. 
     ! If the photon crosses a horizontal boundary of the model, inject 
     ! it into the opposite side. 
     ! If the photon reaches TOA, it exits the simulation.  
     ! If the photon reaches the surface, it is either absorbed or 
     ! reflected.
     ! +-----------------------------------------------------------------
     ! ******************************************************************

     !-------------------------------------------------------------------
     ! Advance photon
     !-------------------------------------------------------------------     
     IF (chatty) PRINT *,''
     IF (chatty) PRINT *,'>>> advance photon >>>'
     IF (chatty) PRINT *,''

     i_old = i
     j_old = j
     k_old = k

     CALL advance_internal(nx, ny, nz, i, j, k, x, y, z, &
             &        i_sign, j_sign, k_sign, a_abs, b_abs, c_abs,&
             &        a_inv, b_inv, c_inv, xunit, yunit, zunit, dist)
   
     IF (chatty) PRINT *,'i=',i,'j=',j,'k=',k
     IF (chatty) PRINT *,'x=',x,'y=',y,'z=',z
     IF (chatty) PRINT *,'a=',a,'b=',b,'c=',c

     !-------------------------------------------------------------------
     ! Save the light path
     !------------------------------------------------------------------- 
     IF (save_path) THEN
        path_layers(k_old,:) = path_layers(k_old,:) + dist
        !IF (k_old .LE. nz_box) THEN
        !   path_boxes(i_old,j_old,k_old,:) = path_boxes(i_old,j_old,k_old,:) + dist
        !ENDIF
     ENDIF

     !-------------------------------------------------------------------
     ! Compute what is left of -ln(RandomNumber).
     !-------------------------------------------------------------------  
     remains = remains - dist*beta_value
     
     IF (chatty) PRINT *,'remains=',remains

     !-------------------------------------------------------------------
     ! If the photon is in a new cell, go back to 1020 and check if an
     ! event has occurred.
     !------------------------------------------------------------------- 
     IF (rough_surface) THEN 
        i_surf_elev = i_surface_elev(i,j)
     ELSE
        i_surf_elev = 1
     ENDIF

     IF ((k .LE. nz).AND.(k .GE. i_surf_elev)) THEN
      
        IF (chatty) PRINT *,'Go back to 1020 to check for atmospheric event'
        GOTO 1020  
   
     !-------------------------------------------------------------------
     ! If the photon reaches TOA, saving the light path is not needed
     ! and go to 4000 to get a new photon.
     !-------------------------------------------------------------------                
     ELSE IF (k .GT. nz) THEN

        IF (chatty) PRINT *,'-------------------'
        IF (chatty) PRINT *,'Photon reached TOA'
        IF (chatty) PRINT *,'-------------------'
        IF (chatty) PRINT *,'GOTO 4000 to get a new photon'
        IF (pauses) READ(*,*)

        GOTO 4000                                      

     !-------------------------------------------------------------------
     ! Else, there is a surface interaction.
     !-------------------------------------------------------------------           
     ELSE

        IF (chatty) PRINT *,'-------------------'     
        IF (chatty) PRINT *,'surface interaction'
        IF (chatty) PRINT *,'-------------------'
        IF (pauses) READ(*,*)

        !update scattering counter 
        IF (save_scats) nscat = nscat + 1
        
        !-------------------------------------------------------------------
        ! Store photon direction before surface interaction 
        ! (needed for new x,y and z definitions af continuation path)
        !-------------------------------------------------------------------          
        a_int = a
        b_int = b
        c_int = c        

        !-------------------------------------------------------------------
        ! Index k cannot be smaller than 1, and at the surface, z is at the
        ! end of the vertical dimension of the cell at the surface
        !-------------------------------------------------------------------         
        k = MAX(k,1)  
        z = 0.D0   

        !-----------------------------------------------------------------------
        ! IF Lambertian surface, then weight the outgoing photon by surface albedo
        !----------------------------------------------------------------------- 
        IF (surface(i,j).EQ.-1) THEN
           CALL find_brdf_albedo(wp)
           weight = weight * wp
        ELSE
           PRINT *,'ERROR: surface should be Lambertian'
           READ(*,*)
        ENDIF

        IF (chatty) PRINT *,'Weight the outgoing photon by the surface albedo'
        IF (chatty) PRINT *,'weight = ',weight

        !--------------------------------------------------
        ! Check if the photon has lost all its weight, 
        ! if yes then inject next photon.
        !--------------------------------------------------        
        IF (weight .LT. eps) THEN
           IF (chatty) PRINT *,'weight < eps'
           IF (chatty) PRINT *,'GOTO 4000 to get a new photon'
           GOTO 4000
        ENDIF

        !-------------------------------------------------------------------
        ! <<<<<<<<<<<<<<<<<<<< Towards the satellite >>>>>>>>>>>>>>>>>>>>>
        !-------------------------------------------------------------------

        !-------------------------------------------------------------------
        ! Get cosines of viewing direction towards the detector for this
        ! location at the surface
        !-------------------------------------------------------------------
        DO ia=1,n_angles    

           !va(ia)=DCOS(DegtoRad*phi_4d(i,j,k,ia))*DSIN(DegtoRad*vartheta_4d(i,j,k,ia))
           !vb(ia)=DSIN(DegtoRad*phi_4d(i,j,k,ia))*DSIN(DegtoRad*vartheta_4d(i,j,k,ia))
           !vc(ia)=DCOS(DegtoRad*vartheta_4d(i,j,k,ia))

           !IF (chatty) PRINT *,'va=',va,'vb=',vb,'vc=',vc

           !-----------------------------------------------------------------------
           ! IF Lambertian surface, then
           !----------------------------------------------------------------------- 
           IF (surface(i,j).EQ.-1) THEN
                  
              !-----------------------------------------------------------------------
              ! Compute the reflected Stokes vector back to the satellite
              !-----------------------------------------------------------------------
              IF (backward_mode) THEN
                 S_back = MATMUL(MATMUL(weight_matrix,As),S0) !Tranpose of As not needed
              ELSE
                 S_back = MATMUL(As,S)
              ENDIF

           ELSE

              !-----------------------------------------------------------------------
              ! Compute the scattering angle and azimuth difference towards detector
              !-----------------------------------------------------------------------
              !phi_view = DegtoRad*phi

              cos_theta_sca = c_int*vc(ia)+DSQRT(MAX(eps,(1.0D+0-c_int**2)))* &
                      & DSQRT(MAX(eps,(1.0D+0-vc(ia)**2)))*DCOS(phi_in_rad-phi_view_rad)

              !-----------------------------------------------------------------------
              ! Get rotation matrices for conversion of polarization reference planes
              !----------------------------------------------------------------------- 
              CALL stokes_rot_mat(L1,L2,cos_theta_sca,c_int,vc(ia),phi_in_rad,&
                         & phi_view_rad,nmat)

              !-----------------------------------------------------------------------
              ! Define the Stokes matrix for surface reflection
              !----------------------------------------------------------------------- 
              !P_surf = 
              PRINT *,'ERROR: Surface should be Lambertian.'
              PRINT *,'surface(i,j)=',surface(i,j)
              READ(*,*)

              !-----------------------------------------------------------------------
              ! Compute the reflected Stokes vector back to the satellite
              ! Forward mode: S = L2 * F_sca * L1 * S
              ! Backward mode: M = M * (L1 * F_sca * L2)^T with L2=L1 and L1=L2 
              ! see Eq. 3.37 of the book of Hovenier et al. (2004).    
              !-----------------------------------------------------------------------
              CALL matmul_stokes_rot(nmat,P_surf,L1,L2,P)

              IF (backward_mode) THEN
                 S_back = MATMUL(MATMUL(weight_matrix,TRANSPOSE(P)),S0)
              ELSE
                 S_back = MATMUL(P,S)     
              ENDIF

           ENDIF

           !-----------------------------------------------------------------------
           ! Compute the Stokes vector contribution
           ! Note that at TOA we divide again by vc 
           !----------------------------------------------------------------------- 
           IF (vc(ia) .LE. 0) THEN
              S_prob(ia,:) = 0.D0 !Signal of surface reflection only travels upwards
           ELSE
              S_prob(ia,:) = S_back * vc(ia) * weight

           ENDIF
           !IF (chatty) PRINT *,'S_prob = ',S_prob

           IF (.NOT. backward_mode) THEN
              IF (nmat .GT. 1) THEN     
                 IF (DABS(S_prob(ia,2)) .GT. eps) THEN
                    PRINT *,'Surface reflection in forward mode but:'
                    PRINT *,'S_prob(ia,2)=',S_prob(ia,2)
                    STOP
                 ENDIF
              ENDIF
           ENDIF

        ENDDO

        !---------------------------------------------------------------------
        ! Add the surface contribution to the radiance arrays
        !---------------------------------------------------------------------
        IF (chatty) PRINT *,'calc_radiances_pol'

        IF (save_scats) THEN
           CALL calc_radiances_and_scats(nx,ny,nz,nmat,n_angles,&
             &     nx_1, ny_1, nz_1, nmat_1, n_angles_1, &
             &     a_int, b_int, c_int, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
             &     TAU_gas,lowest_k_of_only_gas, &
             &     blocked_sensor_path, precompute_detector_paths,&
             &     va, vb, vc,k_detector, &
             &     S_prob,Stokes,nscat,         &
             &     CONTRIB_STOKES_3D,                &
             &     SINGSCAT_CONTRIB_STOKES_3D,                &
             &     MULTSCAT_CONTRIB_STOKES_3D,                &
             &     Stokes_sing,Stokes_mult, &
             &     k_register,backward_mode,i_inject,j_inject, &
             &     i_store_singscat,j_store_singscat,k_store_singscat)

        ELSE

           IF (save_path) THEN
              CALL calc_radiances_and_paths_pol(nx, ny, nz, nmat, n_angles, &
               &     nx_2, ny_2, nz_2, nmat_2, n_angles_2, &
               &     a_int, b_int, c_int, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
               &     TAU_gas,lowest_k_of_only_gas, &
               &     blocked_sensor_path, precompute_detector_paths,&
               &     va, vb, vc, k_detector, &
               &     S_prob,Stokes,path_layers,       &
               !&     nz_box,weight,path_boxes,              &
               &     weighted_path_store, &!weighted_boxpath_store,
               &     sum_of_weights_store,k_register,backward_mode,i_inject,j_inject)
           ELSE
              CALL calc_radiances_pol(nx, ny, nz, nmat, n_angles, &
                &     a_int, b_int, c_int, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
                &     TAU_gas,lowest_k_of_only_gas, &
                &     blocked_sensor_path,precompute_detector_paths, &
                &     va, vb, vc,k_detector, &
                &     S_prob,Stokes,k_register,backward_mode,i_inject,j_inject)
           ENDIF

        ENDIF

        !-------------------------------------------------------------------
        ! <<<<<<<<<<<<<<<<<<<<<< Continuation path >>>>>>>>>>>>>>>>>>>>>>>
        !-------------------------------------------------------------------        
           
        !-------------------------------------------------------------------
        ! IF Lambertian surface, then
        !------------------------------------------------------------------- 
        IF (surface(i,j).EQ.-1) THEN

        !-------------------------------------------------------------------
        ! Find random photon direction after Lambertian surface reflection 
        ! via the rejection method
        !-------------------------------------------------------------------
           !CALL ground_lamb(a,b,c) 
           !phi_out_rad = DATAN2(b,a)

           CALL ground_lamb2(a,b,c,phi_out_rad) !this gives us identical results for fw and bw mc!

           !-------------------------------------------------------------------
           ! Make sure the photon is pointing upwards
           !-------------------------------------------------------------------
           c = DABS(c)              
           
           !-----------------------------------------------------------------------
           ! Compute the reflected Stokes vector for continuation path
           !-----------------------------------------------------------------------      
           IF (backward_mode) THEN
              weight_matrix = MATMUL(weight_matrix,As)
              S(1) = weight_matrix(1,1)
              IF (nmat .GT. 1) THEN
                 S(2) = weight_matrix(1,2)
                 S(3) = -weight_matrix(1,3) !Note the minus sign (needed due to left multiplicatoni by Delta3) 
                 IF (nmat .GT. 3) THEN
                    S(4) = weight_matrix(1,4)
                 ENDIF
              ENDIF
           ELSE
              S = MATMUL(As,S)
           ENDIF

        ELSE

           !-------------------------------------------------------------------
           ! Find random photon direction after surface reflection 
           !-------------------------------------------------------------------                
           !a =
           !b =
           !c =
           !phi_out =
           
           PRINT *,'ERROR: Surface should be Lambertian.'
           PRINT *,'surface(i,j)=',surface(i,j)
           READ(*,*)

           !-----------------------------------------------------------------------
           ! Compute scattering angle and azimuth difference for continuation path
           !-----------------------------------------------------------------------
           cos_theta_sca = c_int*c+DSQRT(MAX(eps,(1.0D+0-c_int**2)))* &
                      & DSQRT(MAX(eps,(1.0D+0-c**2)))*DCOS(phi_in_rad-phi_out_rad)

           !-----------------------------------------------------------------------
           ! Get rotation matrices for conversion of polarization reference planes
           !----------------------------------------------------------------------- 
           CALL stokes_rot_mat(L1,L2,cos_theta_sca,c_int,c,phi_in_rad,&
                   & phi_out_rad,nmat)
                   
           !-----------------------------------------------------------------------
           ! Define the Stokes matrix for surface reflection
           !----------------------------------------------------------------------- 
           !P_surf = 

           !-----------------------------------------------------------------------  
           ! Compute the reflected Stokes vector for the continuation path  
           ! Forward mode: S = L2 * F_sca * L1 * S
           ! Backward mode: M = M * (L1 * F_sca * L2)^T with L2=L1 and L1=L2 
           ! see Eq. 3.37 of the book of Hovenier et al. (2004).    
           !-----------------------------------------------------------------------
           CALL matmul_stokes_rot(nmat,P_surf,L1,L2,P)

           IF (backward_mode) THEN
              weight_matrix = MATMUL(weight_matrix,TRANSPOSE(P))
              S(1) = weight_matrix(1,1)
              IF (nmat .GT. 1) THEN
                 S(2) = weight_matrix(1,2)
                 S(3) = -weight_matrix(1,3) !Note the minus sign (needed due to left multiplicatoni by Delta3) 
                 IF (nmat .GT. 3) THEN
                    S(4) = weight_matrix(1,4)
                 ENDIF
              ENDIF
           ELSE
              S = MATMUL(P,S)
           ENDIF

        ENDIF

        !-----------------------------------------------------------------------
        ! Normalize reflected Stokes vector for continuation path
        !-----------------------------------------------------------------------         
        IF (S(1) .GT. eps) THEN
          S = S/S(1)
        ENDIF

        IF (nmat .EQ. 4) THEN
           IF (SQRT(S(2)**2+S(3)**2+S(4)**2) .GT. S(1)+stokes_criterium_error) THEN
           PRINT *,'ERROR in surface reflection: SQRT(S(2)**2+S(3)**2+S(4)**2) .GT. S(1)'
           ENDIF 
        ELSE IF (nmat .EQ. 3) THEN
           IF (SQRT(S(2)**2+S(3)**2) .GT. S(1)+stokes_criterium_error) THEN
              PRINT *,'ERROR in surface reflection: SQRT(S(2)**2+S(3)**2) .GT. S(1)'
           ENDIF                 
        ENDIF

        IF (backward_mode) THEN
           weight_matrix = weight_matrix/weight_matrix(1,1)
        ENDIF

        !---------------------------------------------------------------------
        ! If the sign of the direction changed, simply reverse x, y and/or z
        !---------------------------------------------------------------------
        IF (A*A_INT .LT. 0.D0) X = MAX(0.D0, XUNIT - X)
        IF (B*B_INT .LT. 0.D0) Y = MAX(0.D0, YUNIT - Y)
        IF (C*C_INT .LT. 0.D0) Z = MAX(0.D0, ZUNIT(K) - Z)

        IF (chatty) PRINT *,'Photon direction after surface reflection: ','a=',a,'b=',b,'c=',c
        IF (chatty) PRINT *,'Photon position after surface reflection: ','i=',i,'j=',j,'k=',k
        IF (chatty) PRINT *,'Photon position after surface reflection: ','x=',x,'y=',y,'z=',z
        IF (pauses) READ(*,*)

        !---------------------------------------------------------------------
        ! Pass the reflected azimuth as input to the continuation path
        !---------------------------------------------------------------------        
        phi_in_rad = phi_out_rad

        !---------------------------------------------------------------------
        !Advance photon after surface reflection
        !---------------------------------------------------------------------        
        IF (chatty) PRINT *,'Advance photon after surface reflection.'
        GOTO 1010

     ENDIF

     !**********************************************************************
     ! +--------------------------------------------------------------------
     ! absorption / scattering event
     ! +--------------------------------------------------------------------
     ! The photon encounters an event in the atmosphere.  
     ! Determine the nature of the event: gas or one of the aerosol types. 
     ! Weight by the single scattering albedo. 
     ! Compute the contribution of the event to the signal.
     ! Finally determine the direction of the photon's continuation path. 
     ! +-------------------------------------------------------------------

     ! +--------------------------------------------------------------------
     !**********************************************************************

       !-------------------------------------------------------------------
       ! Randomly pick scatt type (index 'is') with probability based on
       ! the relative extinction of the scatt type (see Eq. 18 of Barker
       ! et al. 2003) 
       !------------------------------------------------------------------- 
  3000 CALL RANDOM_NUMBER(rtemp)

       is=1
       work=ext_by_species(1,i,j,k) !normalized beta, between 0 and 1

       DO WHILE((work.LE.rtemp).AND.(is.LT.n_types+1))
          is=is+1
          work=work+ext_by_species(is,i,j,k)
       ENDDO
       
       IF (chatty) PRINT *,'scatt type index "is" =',is

       !------------------------------------------------------------
       ! **********************************************************
       ! >>>>>>>>>>>>>>>>>>>>>> Absorption <<<<<<<<<<<<<<<<<<<<<<<<
       ! **********************************************************
       !------------------------------------------------------------ 

       !********** Gas absorption ******************
       IF (is .EQ. n_types+1) THEN       

          weight = weight * ssalb_gas(k)

       !********** Cloud or aerosol absorption *****
       ELSE IF (is .LT. n_types+1) THEN  

         weight = weight * new_scatt_data_pol(is)%ssalb_aer

       !********************************************     
       ELSE

         PRINT *,'is=',is,',n_types=',n_types
         PRINT *,'ERROR: n_types+2 not used anymore'
         STOP

       ENDIF

       !--------------------------------------------------
       ! Check if the photon has lost all its weight, 
       ! if yes then inject next photon.
       ! If not then GOTO 1010 to advance photon
       !--------------------------------------------------        
       IF (weight .LT. eps) THEN
          IF (chatty) PRINT *,'weight < eps'
          IF (chatty) PRINT *,'GOTO 4000 to get a new photon'
          GOTO 4000
       ENDIF

       !------------------------------------------------------------
       ! **********************************************************
       ! >>>>>>>>>>>>>>>>>>>>>> Scattering <<<<<<<<<<<<<<<<<<<<<<<<
       ! **********************************************************
       !------------------------------------------------------------               
       IF (chatty) PRINT *,'Scattering'     

       !update scattering counter
       IF (save_scats) nscat = nscat + 1

       !********** Rayleigh scattering **********
       IF (is .EQ. n_types+1) THEN
         
          IF (chatty) PRINT *,'Rayleigh scattering'
           
       !************ Mie scattering *************  
       ELSE IF (is .LT. n_types+1) THEN
         
         IF (chatty) PRINT *,'Mie scattering'

       ENDIF
       !*****************************************

       IF (pauses) READ(*,*)

       !-------------------------------------------------------------------
       ! <<<<<<<<<<<<<<<<<<<< Towards the satellite >>>>>>>>>>>>>>>>>>>>>
       !-------------------------------------------------------------------
       IF (chatty) PRINT *,'<<<<< Towards the satellite >>>>>'

       !-----------------------------------------------------------------------
       ! Start loop over viewing zenith angles towards the satellite
       !-----------------------------------------------------------------------
       DO ia=1,n_angles   

         !-----------------------------------------------------------------------
         ! Compute the scattering angle and azimuth difference towards satellite
         !-----------------------------------------------------------------------
         !phi_view = DegtoRad*phi

         cos_theta_sca = c*vc(ia)+DSQRT(MAX(eps,(1.0D+0-c**2)))* &
                & DSQRT(MAX(eps,(1.0D+0-vc(ia)**2)))*DCOS(phi_in_rad-phi_view_rad)

         !-----------------------------------------------------------------------
         ! Get rotation matrices for conversion of polarization reference planes
         !----------------------------------------------------------------------- 
         CALL stokes_rot_mat(L1,L2,cos_theta_sca,c,vc(ia),phi_in_rad,&
                    & phi_view_rad,nmat)

         !-------------------------------------------------------------------
         ! Calculate scattering matrix towards satellite
         !-------------------------------------------------------------------

         !********** Rayleigh scattering **********
         IF (is .EQ. n_types+1) THEN
  
           CALL ray_scat_matrix(nmat,cos_theta_sca,depol,F_sca)

         !************ Mie scattering *************  
         ELSE

           CALL mie_scat_matrix(nmat,cos_theta_sca,new_scatt_data_pol(is),new_scatt_data_pol(is)%n_scat_angles,F_sca)   
         
         ENDIF
         !*****************************************

         !-----------------------------------------------------------------------
         ! Compute the reflected Stokes vector back to the satellite
         ! Forward mode: S = L2 * F_sca * L1 * S
         ! Backward mode: M = M * (L1 * F_sca * L2)^T with L2=L1 and L1=L2 
         ! see Eq. 3.37 of the book of Hovenier et al. (2004).            
         !----------------------------------------------------------------------- 
         CALL matmul_stokes_rot(nmat,F_sca,L1,L2,P)         

         IF (backward_mode) THEN
            S_back = MATMUL(MATMUL(weight_matrix,TRANSPOSE(P)),S0)
         ELSE
            S_back = MATMUL(P,S)     
         ENDIF

         !-----------------------------------------------------------------------
         ! Compute the Stokes vector contribution
         !----------------------------------------------------------------------- 

         !Stokes vector
         S_prob(ia,:) = S_back * weight 

         !Neglect extremely small contributions
         IF (S_prob(ia,1) .LT. 1.0D-100) THEN
            S_prob(ia,:)=0.D0
         ENDIF

         !Safety check
         IF (nmat .EQ. 4) THEN
            IF (SQRT(S_prob(ia,2)**2+S_prob(ia,3)**2+S_prob(ia,4)**2) .GT. S_prob(ia,1)+stokes_criterium_error) THEN
               PRINT *,'ERROR in atmospheric event, to satellite:'
               PRINT *,'(SQRT(S_prob(ia,2)**2+S_prob(ia,3)**2+S_prob(ia,4)**2)='
               PRINT *,SQRT(S_prob(ia,2)**2+S_prob(ia,3)**2+S_prob(ia,4)**2)
               PRINT *,'(SQRT(S_prob(ia,2)**2+S_prob(ia,3)**2+S_prob(ia,4)**2) .GT. S_prob(ia,1))'
               PRINT *,'ia=',ia
               PRINT *,'S_prob(ia,:) =',S_prob(ia,:)
               READ(*,*)
               PRINT *,'F_sca=',F_sca
               PRINT *,'L1=',L1
               PRINT *,'L2=',L2
               PRINT *,'is=',is
               PRINT *,'cos_theta_sca=',cos_theta_sca
               !PRINT *,'nscat=',nscat
               READ(*,*)
            ENDIF
         ELSE IF (nmat .EQ. 3) THEN
            IF (SQRT(S_prob(ia,2)**2+S_prob(ia,3)**2) .GT. S_prob(ia,1)+stokes_criterium_error) THEN
               PRINT *,'ERROR in atmospheric event, to satellite:'
               PRINT *,'(SQRT(S_prob(ia,2)**2+S_prob(ia,3)**2)=',SQRT(S_prob(ia,2)**2+S_prob(ia,3)**2) 
               PRINT *,'(SQRT(S_prob(ia,2)**2+S_prob(ia,3)**2) .GT. S_prob(ia,1))'
               PRINT *,'ia=',ia
               PRINT *,'S_prob(ia,:) =',S_prob(ia,:)
               READ(*,*)
               PRINT *,'F_sca=',F_sca
               PRINT *,'L1=',L1
               PRINT *,'L2=',L2
               PRINT *,'is=',is
               PRINT *,'cos_theta_sca=',cos_theta_sca
               !PRINT *,'nscat=',nscat
               READ(*,*)
            ENDIF                  
         ENDIF

       !-----------------------------------------------------------------------
       ! End loop over viewing zenith angles towards the satellite
       !-----------------------------------------------------------------------
       ENDDO

       !------------------------------------------------------------------------
       ! Calculate the radiance and Stokes contributions for each viewing angle
       !------------------------------------------------------------------------
       IF (chatty) PRINT *,'calc_radiances_pol'

       IF (save_scats) THEN
           CALL calc_radiances_and_scats(nx,ny,nz,nmat,n_angles,&
             &     nx_1, ny_1, nz_1, nmat_1, n_angles_1, &
             &     a, b, c, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
             &     TAU_gas,lowest_k_of_only_gas, &
             &     blocked_sensor_path, precompute_detector_paths,&
             &     va, vb, vc,k_detector, &
             &     S_prob,Stokes,nscat,         &
             &     CONTRIB_STOKES_3D,                &
             &     SINGSCAT_CONTRIB_STOKES_3D,                &
             &     MULTSCAT_CONTRIB_STOKES_3D,                &
             &     Stokes_sing,Stokes_mult,&
             &     k_register,backward_mode,i_inject,j_inject, &
             &     i_store_singscat,j_store_singscat,k_store_singscat)
       ELSE

          IF (save_path) THEN
             CALL calc_radiances_and_paths_pol(nx, ny, nz, nmat, n_angles, &
                &     nx_2, ny_2, nz_2, nmat_2, n_angles_2, &                     
                &     a, b, c, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
                &     TAU_gas,lowest_k_of_only_gas, &
                &     blocked_sensor_path, precompute_detector_paths,&
                &     va, vb, vc,k_detector, &
                &     S_prob,Stokes,path_layers,      &
                !&     nz_box,weight,path_boxes,              &
                &     weighted_path_store,&!weighted_boxpath_store,
                &     sum_of_weights_store,k_register,backward_mode,i_inject,j_inject)
          ELSE
             CALL calc_radiances_pol(nx, ny, nz, nmat, n_angles, &
                &     a, b, c, i, j, k, x, y, z, xunit, yunit, zunit, beta, &
                &     TAU_gas,lowest_k_of_only_gas, &
                &     blocked_sensor_path, precompute_detector_paths,&
                &     va, vb, vc,k_detector, &
                &     S_prob,Stokes,k_register,backward_mode,i_inject,j_inject)
          ENDIF

       ENDIF

       !above subroutines are only used to calculate radiance
       !a,b,c,i,j,k,x,y,z did not change (are reset to input values at the end
       !of the subroutine)

       !-----------------------------------------------------------------------
       ! <<<<<<<<<<<<<<<<<<<<<<<<< Continuation path >>>>>>>>>>>>>>>>>>>>>>>>>>
       !-----------------------------------------------------------------------
       IF (chatty) PRINT *,'<<<<< Continuation path >>>>>'
       IF (pauses) READ(*,*)

       !-----------------------------------------------------------------------
       ! Sample the cosine of the scattering and the rotation angle for the
       ! continuation path.
       !----------------------------------------------------------------------- 
       cos_theta_in = c

       IF (is == n_types+1) THEN

         !CALL rayl_scatter_pol(nmat,depol,S,cos_theta_in,cos_theta_out,dphi_rad,cos_theta_sca)     
         !CALL rayl_scatter_pol(backward_mode,nmat,depol,S0,cos_theta_in,cos_theta_out,dphi_rad,cos_theta_sca) !No pol taken into account
         CALL rayl_scatter_pol(backward_mode,nmat,depol,S,cos_theta_in,cos_theta_out,dphi_rad,cos_theta_sca)

       ELSE
         
         PDFerror = .false.
         
         !CALL mie_scatter_pol_bound(nmat,PDFerror,new_scatt_data_pol(is),new_scatt_data_pol(is)%n_scat_angles,&
         !                            & S,cos_theta_in,cos_theta_out,dphi_rad,cos_theta_sca)         
         !CALL mie_scatter_pol_bound(backward_mode,nmat,PDFerror,new_scatt_data_pol(is),new_scatt_data_pol(is)%n_scat_angles,&
         !                           & S0,cos_theta_in,cos_theta_out,dphi_rad,cos_theta_sca) !No pol taken into account
         CALL mie_scatter_pol_bound(backward_mode,nmat,PDFerror,new_scatt_data_pol(is),new_scatt_data_pol(is)%n_scat_angles,&
                                    & S,cos_theta_in,cos_theta_out,dphi_rad,cos_theta_sca)


         IF (PDFerror) THEN
            PRINT*,'S=',S
            PRINT*,'weight=',weight
            PRINT*,'i_photon=',i_photon
         ENDIF
       ENDIF        

       !-----------------------------------------------------------------------
       ! Get rotation matrices for conversion of polarization reference planes
       !----------------------------------------------------------------------- 
       phi_out_rad = phi_in_rad + dphi_rad

       CALL stokes_rot_mat(L1,L2,cos_theta_sca,cos_theta_in,cos_theta_out,phi_in_rad,&
               & phi_out_rad,nmat)

       !-----------------------------------------------------------------------
       ! Get the Stokes matrix for scattering for continuation path
       !----------------------------------------------------------------------- 
       IF (is .EQ. n_types+1) THEN
          
          CALL ray_scat_matrix(nmat,cos_theta_sca,depol,F_sca)
       
       ELSE

          CALL mie_scat_matrix(nmat,cos_theta_sca,new_scatt_data_pol(is),new_scatt_data_pol(is)%n_scat_angles,F_sca)

       ENDIF

       !-----------------------------------------------------------------------
       ! Compute the scattered Stokes vector for continuation path
       ! Forward mode: S = L2 * F_sca * L1 * S
       ! Backward mode: M = M * (L1 * F_sca * L2)^T with L2=L1 and L1=L2 
       ! see Eq. 3.37 of the book of Hovenier et al. (2004).                 
       !----------------------------------------------------------------------- 
       CALL matmul_stokes_rot(nmat,F_sca,L1,L2,P)

       IF (chatty) PRINT *,'S=',S

       IF (backward_mode) THEN
          weight_matrix = MATMUL(weight_matrix,TRANSPOSE(P))
          S(1) = weight_matrix(1,1)
          IF (nmat .GT. 1) THEN
             S(2) = weight_matrix(1,2)
             S(3) = -weight_matrix(1,3) !Note the minus sign (needed due to left multiplicatoni by Delta3) 
             IF (nmat .GT. 3) THEN
                S(4) = weight_matrix(1,4)
             ENDIF
          ENDIF
       ELSE                  
          S = MATMUL(P,S)
       ENDIF       

       !Correction if polarization is not taken into account in theta sampling above (S0 used)
       !The line below is valid for both forward and backward mode
       !weight = weight * S(1) / P(1,1)

       !-----------------------------------------------------------------------
       ! Normalize scattered Stokes vector for continuation path
       !-----------------------------------------------------------------------         
       S = S/S(1)

       IF (chatty) PRINT *,'S normalized =',S

       !Neglect extremely small contributions
       IF (S(1) .LT. 1.0D-100) THEN
          S(:)=0.D0
       ENDIF

       !Safety check
       IF (nmat .EQ. 4) THEN
          IF (SQRT(S(2)**2+S(3)**2+S(4)**2) .GT. S(1)+stokes_criterium_error) THEN
             PRINT *,'ERROR in atmospheric event, continuation path: SQRT(S(2)**2+S(3)**2+S(4)**2) .GT. S(1)'
             PRINT *,'SQRT(S(2)**2+S(3)**2+S(4)**2)=',SQRT(S(2)**2+S(3)**2+S(4)**2)
             PRINT *,'S normalized =',S
             READ(*,*)
             PRINT *,'F_sca=',F_sca 
             PRINT *,'L1=',L1
             PRINT *,'L2=',L2
             PRINT *,'is=',is
             PRINT *,'cos_theta_sca=',cos_theta_sca
             PRINT *,'cos_theta_in',cos_theta_in
             PRINT *,'cos_theta_out',cos_theta_out
             PRINT *,'phi_in_rad',phi_in_rad
             PRINT *,'phi_out_rad',phi_out_rad
             !PRINT *,'nscat=',nscat
             READ(*,*)
          ENDIF
       ELSE IF (nmat .EQ. 3) THEN
          IF (SQRT(S(2)**2+S(3)**2) .GT. S(1)+stokes_criterium_error) THEN
             PRINT *,'ERROR in atmospheric event, continuation path: SQRT(S(2)**2+S(3)**2) .GT. S(1)'
             PRINT *,'SQRT(S(2)**2+S(3)**2)=',SQRT(S(2)**2+S(3)**2)
             PRINT *,'S normalized =',S
             READ(*,*)
             PRINT *,'F_sca=',F_sca
             PRINT *,'L1=',L1
             PRINT *,'L2=',L2
             PRINT *,'is=',is
             PRINT *,'cos_theta_sca=',cos_theta_sca
             PRINT *,'cos_theta_in',cos_theta_in
             PRINT *,'cos_theta_out',cos_theta_out
             PRINT *,'phi_in_rad',phi_in_rad
             PRINT *,'phi_out_rad',phi_out_rad
             !PRINT *,'nscat=',nscat
             READ(*,*)
          ENDIF
       ENDIF

        IF (backward_mode) THEN
           weight_matrix = weight_matrix/weight_matrix(1,1)
        ENDIF

       !---------------------------------------------------------------------- 
       ! Update the direction of the photon according to the scattering angle
       !----------------------------------------------------------------------
       IF (chatty) PRINT *,'Update photon direction.'

       A_BAK = A
       B_BAK = B
       C_BAK = C          

       a=DCOS(phi_out_rad)*SQRT(MAX(0.D0,1.D0-cos_theta_out**2))
       b=DSIN(phi_out_rad)*SQRT(MAX(0.D0,1.D0-cos_theta_out**2))
       c=cos_theta_out !positive upwards                      
     
       !---------------------------------------------------------------------
       ! If the sign of the direction changed, simply reverse x, y and/or z
       !---------------------------------------------------------------------
       IF (A*A_BAK .LT. 0.D0) X = MAX(0.D0, XUNIT - X)
       IF (B*B_BAK .LT. 0.D0) Y = MAX(0.D0, YUNIT - Y)
       IF (C*C_BAK .LT. 0.D0) Z = MAX(0.D0, ZUNIT(K) - Z)   

       IF (chatty) PRINT *,'i=',i,'j=',j,'k=',k
       IF (chatty) PRINT *,'x=',x,'y=',y,'z=',z
       IF (chatty) PRINT *,'a=',a,'b=',b,'c=',c

       !---------------------------------------------------------------------
       ! Pass the reflected azimuth as input to the continuation path
       !---------------------------------------------------------------------              
       phi_in_rad = phi_out_rad

       !---------------------------------------------------------------------
       ! Advance photon after scattering 
       ! GOTO 1000 for new photon path distance and new direction
       !---------------------------------------------------------------------   
       IF (chatty) PRINT *,'Advance photon after scattering.'       
       GOTO 1000    

  !-------------------------------------------------------------------
  ! Photon exits simulation
  !-------------------------------------------------------------------  
  4000 ENDDO 

  !*******************************************************************
  ! +-----------------------------------------------------------------
  ! End of simulation.
  ! +-----------------------------------------------------------------
  ! All photons have been put through the model. 
  ! Scale radiances and exit the subroutine.
  ! +-----------------------------------------------------------------
  !*******************************************************************
  IF (chatty) PRINT *,''
  IF (chatty) PRINT *,'End of simulation.'
  IF (chatty) PRINT *,'All photons have been put through the model.'
  IF (chatty) PRINT *,'Scale radiances and exit the subroutine.'

  !----------------------------------------------------------------------
  ! In forward mode, multiply by the mu_0/mu factor, originating from
  ! the conversion between the horizontal source and detector areas
  ! associated with the same scattering voume element
  !----------------------------------------------------------------------
  IF (.NOT. backward_mode) THEN

     !loop over zenith angles
     DO ia=1,n_angles

        !compute mu0/mu factor
        mu0mu_factor = DABS(DCOS(vartheta0_rad))/DABS(DCOS(vartheta_rad(ia)))

        !multiply by mu0/mu
        Stokes(:,:,:,ia) = Stokes(:,:,:,ia) * mu0mu_factor
        IF (save_scats) THEN
           Stokes_sing(:,:,:,ia) = Stokes_sing(:,:,:,ia) * mu0mu_factor
           Stokes_mult(:,:,:,ia) = Stokes_mult(:,:,:,ia) * mu0mu_factor
           CONTRIB_STOKES_3D(:,:,:,:,ia) = CONTRIB_STOKES_3D(:,:,:,:,ia) * mu0mu_factor
           SINGSCAT_CONTRIB_STOKES_3D(:,:,:,:,ia) = SINGSCAT_CONTRIB_STOKES_3D(:,:,:,:,ia) * mu0mu_factor
           MULTSCAT_CONTRIB_STOKES_3D(:,:,:,:,ia) = MULTSCAT_CONTRIB_STOKES_3D(:,:,:,:,ia) * mu0mu_factor
        ENDIF
     ENDDO
  ENDIF

  RETURN

CONTAINS

!---------------------------------------------------------------------
!********************** INTERNAL SUBROUTINES *************************
!---------------------------------------------------------------------

! NOTE THAT VARIABLES ARE SHARED WITH PROGRAM ABOVE!
! That is, passing input and output arguments explicitly is not 
! required, but it is good practice for code readability.

!+----------------------------------------------------------------
! SUBROUTINE find_brdf_albedo
!+----------------------------------------------------------------
! This subroutine returns the surface albedo at location (i,j)
!+----------------------------------------------------------------
  SUBROUTINE find_brdf_albedo(wp)

    ! Output ------------------------------------------------------
    REAL*8,intent(out)   :: wp !Variable to store surface albedo (-)
    
    !--------------------------------------------------------------
    ! Get surface class
    !--------------------------------------------------------------
    !surface_class_brdf=surface(i,j)
        
    !--------------------------------------------------------------
    ! IF the surface is Lambertian, put surface albedo in WP
    !--------------------------------------------------------------
    IF (surface(i,j)==-1) THEN
       
       wp = surface_alb(i,j)
       
    !-------------------------------------------------------------------
    ! ELSE stop the code, this version only handles Lambertian surfaces
    !-------------------------------------------------------------------   
    ELSE
       PRINT *,'ERROR: code can only handle Lambertian surfaces.'
       STOP
    ENDIF   

    RETURN

    END SUBROUTINE find_brdf_albedo


!+----------------------------------------------------------------
! SUBROUTINE ADVANCE_internal
!+----------------------------------------------------------------
! Advance the photon to the next cell and calculate DIST, the 
! distance it travelled to get there.  If the photon crosses a
! horizontal boundary of the model, inject it into the opposite
! side.
!+----------------------------------------------------------------  
  SUBROUTINE ADVANCE_internal(NX, NY, NZ, I, J, K, X, Y, Z, &
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
END SUBROUTINE ADVANCE_internal

END SUBROUTINE monte_carlo
