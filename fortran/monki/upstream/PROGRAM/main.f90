! SPDX-License-Identifier: MIT
!
!******************************************************************************
!*
!* MONKI (Monte Carlo KNMI) -- Cartesian Version
!*
!* Software version: 1.0.1
!* Release date:     2026-09-13
!* Release status:   Final publication release
!*
!* AUTHORSHIP AND DEVELOPMENT
!* --------------------------
!* Author, creator, and lead developer:
!*   Victor J.H. Trees
!*
!* Developed while employed by TU Delft and working at KNMI.
!*
!* Contact:
!*   victor.trees@knmi.nl
!*   v.j.h.trees@tudelft.nl
!*
!* Copyright (c) 2026 Victor J.H. Trees
!*
!* IMPORTANT: SCIENTIFIC CITATION AND CONTACT
!* ------------------------------------------
!* MONKI is scientific research software developed by Victor J.H. Trees.
!* Proper acknowledgement of the software and its lead developer is
!* important for recognising the scientific and software-development work
!* underlying MONKI.
!*
!* If MONKI contributes to research, analysis, figures, results, software,
!* reports, presentations, or publications, users are strongly requested
!* to cite both:
!*
!*   1. Trees et al. (2026), "MONKI: a three-dimensional Monte Carlo
!*      simulator for transfer of polarised light in planetary
!*      atmospheres", Geoscientific Model Development.
!*
!*   2. The DOI of the specific MONKI software release used.
!*
!* Please do not refer to MONKI as an unidentified or generic Monte Carlo
!* model. Clearly identify MONKI by name and acknowledge Victor J.H. Trees
!* as its creator and lead developer.
!*
!* Users intending to apply MONKI in a substantial scientific study,
!* modify or extend the model, incorporate it into another model, or use
!* it for operational, mission-related, or commercial work are strongly
!* encouraged to contact Victor J.H. Trees at an early stage.
!*
!* Contact supports correct use and interpretation of MONKI, avoids
!* duplication of development, and provides an opportunity to coordinate
!* improvements and discuss scientific collaboration. Collaboration is
!* warmly welcomed.
!*
!* These requests reflect accepted scientific practice and are separate
!* from the legal conditions of the MIT License. Citation, contact, and
!* collaboration are not mandatory conditions of the MIT License.
!*
!* DESCRIPTION
!* -----------
!* Three-dimensional radiative transfer simulation of sunlight reflected
!* by an atmosphere-surface system using the Monte Carlo technique. MONKI
!* fully accounts for the polarisation state of light (Stokes parameters
!* I, Q, U, and V) at all scattering orders.
!*
!* Atmosphere:
!*   - Rayleigh scattering by gaseous molecules
!*   - Mie scattering by spherical droplets
!*   - Absorption by gases
!*   - Absorption by spherical droplets
!*
!* Surface:
!*   - Lambertian, isotropically reflecting and fully depolarising,
!*     with a specified surface albedo
!*
!* Domain:
!*   - Cartesian grid
!*   - Cyclic boundaries in the horizontal x- and y-directions
!*   - nx grid cells in the x-direction
!*   - ny grid cells in the y-direction
!*   - nz+1 grid levels in the z-direction
!*
!* LICENSE
!* -------
!* MONKI is distributed under the MIT License. Permission is granted to
!* use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!* copies of the software, subject to the conditions of the MIT License.
!*
!* The above copyright notice and the MIT License notice must be included
!* in all copies or substantial portions of the software.
!*
!* MONKI is provided "as is", without warranty of any kind, express or
!* implied, including but not limited to the warranties of merchantability,
!* fitness for a particular purpose, and non-infringement.
!*
!* See the LICENSE file in the root directory for the complete MIT License.
!*
!* THIRD-PARTY COMPONENTS
!* ----------------------
!* The ncdf_utilities.f90 included with MONKI originates from ECSIM and is
!* distributed under the GNU Lesser General Public License. Its original
!* copyright and license notices are retained in the source file.
!* Third-party components remain subject to their respective licenses.
!*
!******************************************************************************

  PROGRAM main
 
  !---------------------------------------------------------------------
  ! Import modules
  !---------------------------------------------------------------------
  USE main_routines_general 
  USE main_routines_cartesian
  USE new_data_types
  USE write_messages
  USE physical_parameters
  USE netcdf
  USE ncdf_utilities

  !---------------------------------------------------------------------
  ! Declare variables
  !---------------------------------------------------------------------
  IMPLICIT NONE

  !Include max size of grid variable arrays
  INCLUDE 'max_incl'

  ! Command line input ---------------------------------------------------------------------
  INTEGER :: n_angles          !Number of viewing zenith angles (-)
  INTEGER :: backward_mc       !Logical to indicate whether code runs in backward mode (-)

  REAL*8 :: asurf              !Surface albedo (-)
  REAL*8 :: theta0             !Solar zenith angle (deg)
  REAL*8 :: phi0               !Solar azimuth angle (deg)
  REAL*8 :: phi                !Viewing azimuth angle (deg)
  
  REAL*8,DIMENSION(:),ALLOCATABLE :: vartheta !Viewing zenith angles (deg)
  
  CHARACTER(len=500) :: infile      !Input filename
  CHARACTER(len=500) :: outfile_ext !Output filename extension
  
  ! Written output ------------------------------------------------------------------------
  INTEGER :: N_ROUNDS             !Number of simulation rounds used (-)
  
  INTEGER*8 :: N_PHOTONS_PER_COLUMN !Number of photons used per column and per simulation round (-)

  REAL*8 :: x_unit                !Dimensions of the grid cells in the x direction (km)
  REAL*8 :: y_unit                !Dimensions of the grid cells in the y direction (km)
  REAL*8 :: wav                   !Wavelength of the light (nm)

  REAL*8,DIMENSION(:),ALLOCATABLE :: z_unit !Dimensions of the grid cells in the z direction (km)
  REAL*8,DIMENSION(:),ALLOCATABLE :: pres   !Pressure at the levels (Pa)
  REAL*8,DIMENSION(:),ALLOCATABLE :: temp   !Temperature at the levels (K)
  REAL*8,DIMENSION(:),ALLOCATABLE :: alt    !Altitude at the levels (m)
  REAL*8,DIMENSION(:),ALLOCATABLE :: bmabs  !Gas absorption optical thickness in the layers (-) 
  REAL*8,DIMENSION(:),ALLOCATABLE :: bmsca  !Gas scattering optical thickness in the layers (-)
  !REAL*8,DIMENSION(:),ALLOCATABLE :: S1_ROUNDS_save !Convergence history of the mean radiance per round (-)
  
  REAL*8,DIMENSION(:,:),ALLOCATABLE :: STOKES_ROUNDS_save !Convergence history of the mean radiance per round (-)
  REAL*8,DIMENSION(:,:),ALLOCATABLE :: vmrgas !Volume mixing ratios of the gases at the levels (ppm)

  REAL*8,DIMENSION(:,:,:),ALLOCATABLE :: DOP !Degree of polarization (-)

  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: baer !Aerosol total optical thicknes (baerabs+baersca) in the cells (-)  
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: SUM_STOKES !Radiance Stokes vector (-)
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: STD_STOKES !Standard deviation or radiance (-)
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: SUM_STOKES_SING !Single scattering contribution to the signal (-)
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: SUM_STOKES_MULT !Multiple scattering contribution to the signal (-)
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: AMF_layer !Layer air mass factor

  REAL*8,DIMENSION(:,:,:,:,:),ALLOCATABLE :: SUM_CONTRIB_STOKES_3D !3D contribution to the signal (-)
  REAL*8,DIMENSION(:,:,:,:,:),ALLOCATABLE :: SUM_SINGSCAT_CONTRIB_STOKES_3D !3D single scattering contribution to the signal (-)
  REAL*8,DIMENSION(:,:,:,:,:),ALLOCATABLE :: SUM_MULTSCAT_CONTRIB_STOKES_3D !3D multiple scattering contribution to the signal (-)

  ! Internal ------------------------------------------------------------------------------
  INTEGER :: nmat,nt,nx,ny,nz,nbands,iband,nscatt_types,nwaves
  INTEGER :: nmat_1,nx_1,ny_1,nz_1,n_angles_1,nmat_2,nx_2,ny_2,nz_2,n_angles_2
  INTEGER :: n_scat_angles,outputnum,nclouds,gas_col_type
  INTEGER :: iangles,is,ia,seed,status,i,j,k,isc,seed_size,seedidx
  INTEGER :: k_detector,iround,k_register
  INTEGER :: lowest_k_of_only_gas,k0,ngases,nlev
  INTEGER :: iw,detector_location,k_source,register_location!,nz_box,boxamfID
  INTEGER :: source_location,detector_location_save
  INTEGER :: stokesID,dopID,amfID,stokessingID,stokesmultID
  INTEGER :: contrib_stokesID,singscat_contrib_stokesID,stdstokesID
  INTEGER :: multscat_contrib_stokesID,bmabsID,bmscaID,baerID,wavID
  INTEGER :: start_clock, end_clock, clock_rate
  INTEGER :: total_seconds, days, hours, minutes, seconds

  INTEGER*8 :: N_PHOTONS_PER_ROUND,N_PHOTONS,n_columns

  INTEGER,PARAMETER :: ibmsca=11,ibmabs=12

  INTEGER,SAVE :: myid  

  !$ INTEGER,EXTERNAL :: omp_get_thread_num               
  !$ INTEGER,EXTERNAL :: omp_get_num_threads

  INTEGER, ALLOCATABLE  :: seedarray(:)

  INTEGER,DIMENSION(:),ALLOCATABLE :: aerosol_present

  INTEGER,DIMENSION(:,:),ALLOCATABLE :: surface,i_surface_elev
  
  INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE :: blocked_sensor_path

  REAL   :: elapsed_seconds,milliseconds

  REAL*8 :: depol_fixed,vartheta0
  REAL*8 :: wave,STD_S1_MAX,tol_error
  REAL*8 :: depol,x_dim,y_dim
  REAL*8 :: wave_min,wave_max,dwave,wav_read,wavnr
  REAL*8 :: STD_S1_MEAN

  REAL*8,PARAMETER :: eps = 1.0D-16

  REAL*8,DIMENSION(:),ALLOCATABLE :: ssalb_gas
  REAL*8,DIMENSION(:),ALLOCATABLE :: BETA_gas

  REAL*8,DIMENSION(:,:),ALLOCATABLE :: surface_alb_p 
  REAL*8,DIMENSION(:,:),ALLOCATABLE :: TAU_gas

  REAL*8,DIMENSION(:,:,:),ALLOCATABLE :: BETA
  REAL*8,DIMENSION(:,:,:),ALLOCATABLE :: SUM_sum_of_weights_store

  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: cloud_abundance
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: ext_by_scatt_species
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: norm_ext_by_scatt_species
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: SUM_STOKES_ROUNDS,SUM_STOKES_ROUNDS_SQ
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: ERROR_STOKES
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: STOKES,STOKES_SING,STOKES_MULT
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: SUM_weighted_path_store
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: mean_path_layer
  
  REAL*8,DIMENSION(:,:,:,:,:),ALLOCATABLE :: CONTRIB_STOKES_3D
  REAL*8,DIMENSION(:,:,:,:,:),ALLOCATABLE :: SINGSCAT_CONTRIB_STOKES_3D
  REAL*8,DIMENSION(:,:,:,:,:),ALLOCATABLE :: MULTSCAT_CONTRIB_STOKES_3D

  !REAL*8,DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: SUM_weighted_boxpath_store
  !REAL*8,DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: mean_path_box,AMF_box

  REAL*8,DIMENSION(max_nmat,max_nx,max_ny,max_na) :: STOKES_p
  
  REAL*8,DIMENSION(max_nmat_1,max_nx_1,max_ny_1,max_na_1) :: STOKES_SING_p,STOKES_MULT_p
  REAL*8,DIMENSION(max_nmat_1,max_nx_1,max_ny_1,0:max_nz_1,max_na_1) :: CONTRIB_STOKES_3D_p
  REAL*8,DIMENSION(max_nmat_1,max_nx_1,max_ny_1,0:max_nz_1,max_na_1) :: SINGSCAT_CONTRIB_STOKES_3D_p
  REAL*8,DIMENSION(max_nmat_1,max_nx_1,max_ny_1,0:max_nz_1,max_na_1) :: MULTSCAT_CONTRIB_STOKES_3D_p 

  REAL*8,DIMENSION(max_nx_2,max_ny_2,max_nz_2,max_na_2) :: weighted_path_store
  !REAL*8,DIMENSION(max_nx_2,max_ny_2,max_nx_2,max_ny_2,max_nz_box_2,max_na_2) :: weighted_boxpath_store
  REAL*8,DIMENSION(max_nx_2,max_ny_2,max_na_2) :: sum_of_weights_store

  CHARACTER(len=2)  :: wave_type
  CHARACTER(len=3)  :: sasurf,stheta0,sphi0,sphi  
  CHARACTER(len=4)  :: cloud_input_type
  CHARACTER(len=5)  :: background
  CHARACTER(len=20)  :: wrong_dim

  CHARACTER(len=500) :: bmsca_file,bmabs_file

  CHARACTER(len=500) :: arg  

  CHARACTER(len=120) :: string,error_str

  CHARACTER(len=500) :: output_nc_file,name
  CHARACTER(len=300) :: monki_home,signal_name,loc_name,loc_with_signal_name

  CHARACTER*30,DIMENSION(:),ALLOCATABLE :: coefile,coelut
 
  LOGICAL :: finished_rad,save_path,chatty,save_scats
  LOGICAL :: read_bmsca_from_file,read_bmabs_from_file,stop_converged
  LOGICAL :: basicinfo,fix_depol,backward_mode,save_convergence_history
  LOGICAL :: precompute_detector_paths

  TYPE(new_scatterer_info_pol),DIMENSION(:),ALLOCATABLE :: new_scatt_info_pol

  !---------------------------------------------------------------------
  !*************************** START SCRIPT ****************************
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Get the total number of threads   
  !---------------------------------------------------------------------
  nt = 1
  !$OMP PARALLEL
  nt = omp_get_num_threads()
  !$OMP END PARALLEL

  !---------------------------------------------------------------------
  ! Determine the required size for the seed array
  !---------------------------------------------------------------------
  CALL RANDOM_SEED(SIZE=seed_size)

  !---------------------------------------------------------------------
  ! Parallel part for setting a unique seed for each thread
  !---------------------------------------------------------------------
  seed = 23274

  !$OMP PARALLEL PRIVATE(myid, seedarray, seedidx)
      ALLOCATE(seedarray(seed_size))     

      myid = omp_get_thread_num()
      seedarray = [(seed + myid + seedidx, seedidx = 1, seed_size)]
      CALL RANDOM_SEED(PUT=seedarray)

      DEALLOCATE(seedarray)        
  !$OMP END PARALLEL

  !---------------------------------------------------------------------
  ! Decide whether the script writes debugging output to screen for this
  ! main file
  !---------------------------------------------------------------------
  chatty = .false.
  basicinfo = .true.
  
  !---------------------------------------------------------------------
  ! Decide whether to save the convergence history of the result
  !---------------------------------------------------------------------  
  save_convergence_history = .true.

  !---------------------------------------------------------------------
  ! Get the start time
  !---------------------------------------------------------------------
  CALL system_clock(start_clock, clock_rate)
  
  !---------------------------------------------------------------------
  ! Get the MONKI home environment variable
  !---------------------------------------------------------------------
  CALL getenv('MONKI_HOME', monki_home)

  !---------------------------------------------------------------------
  ! Get command line arguments
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Get command line arguments'

  !Surface albedo
  CALL get_command_argument(1, arg)
  READ(arg,*) asurf
  IF (chatty) PRINT *,'asurf = ',asurf

  !Solar zenith angle
  CALL get_command_argument(2, arg)
  READ(arg,*) theta0
  IF (chatty) PRINT *,'theta0 = ',theta0

  !Solar azimuth angle
  CALL get_command_argument(3, arg)
  READ(arg,*) phi0
  IF (chatty) PRINT *,'phi0 = ',phi0

  !Viewing azimuth angle
  CALL get_command_argument(4, arg)
  READ(arg,*) phi
  IF (chatty) PRINT *,'phi = ',phi

  !Number of viewing zenith angles
  CALL get_command_argument(5, arg)
  READ(arg,*) n_angles
  IF (chatty) PRINT *,'n_angles = ',n_angles

  !Viewing zenith angles
  ALLOCATE(vartheta(n_angles))
  vartheta=0.0D0
  DO iangles=1,n_angles
     CALL get_command_argument(5+iangles, arg)
     READ(arg,*) vartheta(iangles)
     IF (chatty) PRINT *,'iangles=',iangles,', vartheta(iangles) = ',vartheta(iangles)
  ENDDO

  !Backward mode logical
  CALL get_command_argument(6+n_angles, arg)
  READ(arg,*) backward_mc
  IF (chatty) PRINT *,'backward_mc = ',backward_mc

  !Name of input file
  CALL get_command_argument(7+n_angles, arg)
  infile = TRIM(ADJUSTL(monki_home)) // '/INPUT/' // TRIM(ADJUSTL(arg))
  IF (chatty) PRINT *,'infile = ',infile

  !Name of output file
  CALL get_command_argument(8+n_angles, arg)
  outfile_ext = TRIM(ADJUSTL(arg))
  IF (chatty) PRINT *,'outfile_ext = ',outfile_ext

  !---------------------------------------------------------------------
  ! Check if the code needs to run in backward Monte Carlo mode
  !---------------------------------------------------------------------
  IF (backward_mc .EQ. 1) THEN
     backward_mode = .true.
  ELSE IF (backward_mc .EQ. 0) THEN
     backward_mode = .false.
  ELSE
     PRINT *,'backward_mc=',backward_mc
     PRINT *,'Error in main.f90: backward_mc must be 0 or 1'
     STOP
  ENDIF

  !---------------------------------------------------------------------
  ! Read the grid size of the scene
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'CALL rdsize'

  !Read dimensions from input file
  CALL rdsize(infile,nx,ny,nz,x_dim,y_dim,nmat,nscatt_types,&
  &   n_scat_angles,wave_type,wave_min,wave_max,dwave,ngases,nclouds)

  !fixed inputs (can possibly be removed later)
  !nz_box = 2
  
  !Number of levels between layers
  nlev = nz + 1

  !Number of wavelengths or wavenumbers
  nwaves = INT((wave_max - wave_min) / dwave) + 1

  !---------------------------------------------------------------------
  ! Check if there is a dimension error
  !---------------------------------------------------------------------  
  wrong_dim = 'none'

  IF (nx .GT. max_nx) THEN 
     wrong_dim = 'nx'
  ELSEIF (ny .GT. max_ny) THEN 
     wrong_dim = 'ny'
  ELSEIF (nz .GT. max_nz) THEN
     wrong_dim = 'nz'     
  ELSEIF (nmat .GT. max_nmat) THEN
     wrong_dim = 'nmat'
  ELSEIF (n_angles .GT. max_na) THEN     
     wrong_dim = 'n_angles'
  !ELSEIF (nz_box .GT. max_nz_box) THEN  
  !   wrong_dim = 'nz_box'
  ELSE
     wrong_dim = 'none'
  ENDIF

  IF (wrong_dim .NE. 'none') THEN
     PRINT *,'Dimension error'
     PRINT *,wrong_dim,' is larger than max value set in max_incl'
     STOP
  ENDIF

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Allocate arrays'

  ALLOCATE(surface(nx,ny))
  ALLOCATE(i_surface_elev(nx,ny))
  ALLOCATE(surface_alb_p(nx,ny))
  ALLOCATE(STOKES(nmat,nx,ny,n_angles))
  ALLOCATE(DOP(nx,ny,n_angles))
  ALLOCATE(SUM_STOKES(nmat,nx,ny,n_angles))
  ALLOCATE(SUM_STOKES_ROUNDS(nmat,nx,ny,n_angles))
  ALLOCATE(SUM_STOKES_ROUNDS_SQ(nmat,nx,ny,n_angles))
  ALLOCATE(STD_STOKES(nmat,nx,ny,n_angles))
  ALLOCATE(ERROR_STOKES(nmat,nx,ny,n_angles))
  ALLOCATE(BETA(nx,ny,nz))
  ALLOCATE(ext_by_scatt_species(nscatt_types+1,nx,ny,nz))
  ALLOCATE(norm_ext_by_scatt_species(nscatt_types+1,nx,ny,nz))
  ALLOCATE(BETA_gas(nz))
  ALLOCATE(TAU_gas(nz,n_angles))
  ALLOCATE(blocked_sensor_path(nx,ny,nz,n_angles))
  ALLOCATE(bmabs(nz))
  ALLOCATE(bmsca(nz))
  ALLOCATE(baer(nscatt_types,nx,ny,nz))
  ALLOCATE(cloud_abundance(nscatt_types,nx,ny,nz))
  ALLOCATE(alt(nlev))
  ALLOCATE(pres(nlev))
  ALLOCATE(temp(nlev))
  ALLOCATE(z_unit(nz))
  ALLOCATE(ssalb_gas(nz))
  ALLOCATE(coefile(nscatt_types))
  ALLOCATE(coelut(nscatt_types))
  ALLOCATE(aerosol_present(nz))
  ALLOCATE(vmrgas(ngases,nlev))

  !---------------------------------------------------------------------
  ! Allocate and initialize scattering matrices elements datatype
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Allocate and initialize scattering matrices elements datatype'

  ALLOCATE(new_scatt_info_pol(1:nscatt_types))
  DO is=1,nscatt_types
     CALL allocate_and_initialize_new_scatterer_info_pol(new_scatt_info_pol(is), &
    &                                                    n_scat_angles)
  ENDDO

  !---------------------------------------------------------------------
  ! Read the input file
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'CALL rdinput'

  CALL rdinput(infile,nx,ny,nz,nscatt_types,ngases,nclouds,            &
    &          vmrgas,pres,temp,alt,                                   &
    &          read_bmsca_from_file,bmsca_file,depol_fixed,            &
    &          gas_col_type,fix_depol,read_bmabs_from_file,bmabs_file, &
    &          background,coefile,coelut,                              &
    &          cloud_input_type,cloud_abundance,                       &
    &          N_PHOTONS_PER_COLUMN,N_ROUNDS,                          &
    &          save_path,save_scats,                                   &
    &          stop_converged,tol_error,detector_location,             &
    &          source_location,register_location)

  !---------------------------------------------------------------------
  ! Allocate large arrays only when necessary
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Allocate large arrays only when necessary'

  IF (save_scats) THEN
     nmat_1     = nmat
     nx_1       = nx
     ny_1       = ny
     nz_1       = nz
     n_angles_1 = n_angles
  ELSE
     nmat_1     = 1
     nx_1       = 1
     ny_1       = 1
     nz_1       = 1
     n_angles_1 = 1
  ENDIF

  ALLOCATE(CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1))
  ALLOCATE(SINGSCAT_CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1))
  ALLOCATE(MULTSCAT_CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1))
  ALLOCATE(SUM_CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1))
  ALLOCATE(SUM_SINGSCAT_CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1))
  ALLOCATE(SUM_MULTSCAT_CONTRIB_STOKES_3D(nmat_1,nx_1,ny_1,0:nz_1,n_angles_1))
  ALLOCATE(STOKES_SING(nmat_1,nx_1,ny_1,n_angles_1))
  ALLOCATE(STOKES_MULT(nmat_1,nx_1,ny_1,n_angles_1))
  ALLOCATE(SUM_STOKES_SING(nmat_1,nx_1,ny_1,n_angles_1))
  ALLOCATE(SUM_STOKES_MULT(nmat_1,nx_1,ny_1,n_angles_1))
  
  IF (save_convergence_history) THEN
     !ALLOCATE(S1_ROUNDS_save(n_rounds))
     ALLOCATE(STOKES_ROUNDS_save(nmat,n_rounds))
  ENDIF

  IF (save_path) THEN
     nmat_2     = nmat
     nx_2       = nx
     ny_2       = ny
     nz_2       = nz
     n_angles_2 = n_angles
  ELSE
     nmat_2     = 1
     nx_2       = 1
     ny_2       = 1
     nz_2       = 1
     n_angles_2 = 1
  ENDIF

  ALLOCATE(SUM_weighted_path_store(nx_2,ny_2,nz_2,n_angles_2))
  !ALLOCATE(SUM_weighted_boxpath_store(nx_2,ny_2,nx_2,ny_2,nz_box_2,n_angles_2))
  ALLOCATE(SUM_sum_of_weights_store(nx_2,ny_2,n_angles_2))
  ALLOCATE(mean_path_layer(nx_2,ny_2,nz_2,n_angles_2))
  ALLOCATE(AMF_layer(nx_2,ny_2,nz_2,n_angles_2))
  !ALLOCATE(mean_path_box(nx_2,ny_2,nx_2,ny_2,nz_box_2,n_angles_2))
  !ALLOCATE(AMF_box(nx_2,ny_2,nx_2,ny_2,nz_box_2,n_angles_2))

  !---------------------------------------------------------------------
  ! Check if there is a dimension error
  !--------------------------------------------------------------------- 
  IF (chatty) PRINT *,'Check if there is a dimension error'

  wrong_dim = 'none'

  IF (save_scats) THEN
     IF (max_nx_1 .NE. max_nx) wrong_dim = 'max_nx_1 .NE. max_nx'
     IF (max_ny_1 .NE. max_ny) wrong_dim = 'max_ny_1 .NE. max_ny'
     IF (max_nz_1 .NE. max_nz) wrong_dim = 'max_nz_1 .NE. max_nz'
     IF (max_nmat_1 .NE. max_nmat) wrong_dim = 'max_nmat_1 .NE. max_nmat'
     IF (max_na_1 .NE. max_na) wrong_dim = 'max_na_1 .NE. max_na'
  ELSE
     IF (max_nx_1 .NE. 1) wrong_dim = 'max_nx_1 .NE. 1'
     IF (max_ny_1 .NE. 1) wrong_dim = 'max_ny_1 .NE. 1'
     IF (max_nz_1 .NE. 1) wrong_dim = 'max_nz_1 .NE. 1'
     IF (max_nmat_1 .NE. 1) wrong_dim = 'max_nmat_1 .NE. 1'
     IF (max_na_1 .NE. 1) wrong_dim = 'max_na_1 .NE. 1'
  ENDIF

  IF (wrong_dim .NE. 'none') THEN
     PRINT *,'Dimension error'
     PRINT *,wrong_dim,' , while in the input file you have set save_scats =',save_scats
     PRINT *,'Please adjust max_incl and compile the program again'
     STOP
  ENDIF

  wrong_dim = 'none'

  IF (save_path) THEN
     IF (max_nx_2 .NE. max_nx) wrong_dim = 'max_nx_2 .NE. max_nx'
     IF (max_ny_2 .NE. max_ny) wrong_dim = 'max_ny_2 .NE. max_ny'
     IF (max_nz_2 .NE. max_nz) wrong_dim = 'max_nz_2 .NE. max_nz'
     IF (max_nmat_2 .NE. max_nmat) wrong_dim = 'max_nmat_2 .NE. max_nmat'
     IF (max_na_2 .NE. max_na) wrong_dim = 'max_na_2 .NE. max_na'
  ELSE
     IF (max_nx_2 .NE. 1) wrong_dim = 'max_nx_2 .NE. 1'
     IF (max_ny_2 .NE. 1) wrong_dim = 'max_ny_2 .NE. 1'
     IF (max_nz_2 .NE. 1) wrong_dim = 'max_nz_2 .NE. 1'
     IF (max_nmat_2 .NE. 1) wrong_dim = 'max_nmat_2 .NE. 1'
     IF (max_na_2 .NE. 1) wrong_dim = 'max_na_2 .NE. 1'
  ENDIF

  IF (wrong_dim .NE. 'none') THEN
     PRINT *,'Dimension error'
     PRINT *,wrong_dim,' , while in the input file you have set save_path =',save_path
     PRINT *,'Please adjust max_incl and compile the program again'
     STOP
  ENDIF

  !---------------------------------------------------------------------
  ! Compute cell sizes in km and compute number of columns
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Compute cell sizes in km and compute number of columns'

  x_unit = x_dim/REAL(nx,8)
  y_unit = y_dim/REAL(ny,8)
  DO k=1,nz
     z_unit(k) = (alt(k+1)-alt(k))/1000.D0
  ENDDO

  n_columns = nx*ny

  !---------------------------------------------------------------------
  ! Specify output filename
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Specify output filename'

  outputnum = 0 !Initialize output file existence logical
  
  WRITE(sasurf,'(I3.3)') NINT(asurf*100)  
  WRITE(stheta0,'(I3.3)') NINT(theta0)
  
  IF (backward_mode) THEN
     WRITE(sphi0,'(I3.3)') NINT(MOD(phi+180.D0,360.D0))
     WRITE(sphi,'(I3.3)') NINT(MOD(phi0+180.D0,360.D0))            
     output_nc_file = TRIM(ADJUSTL(outfile_ext))                         &
     &              //'_as_'//sasurf//'_th_'//stheta0                    &
     &              //'_ph0_'//sphi0//'_ph_'//sphi//'.nc'
  ELSE
     WRITE(sphi0,'(I3.3)') NINT(phi0)
     WRITE(sphi,'(I3.3)') NINT(phi)     
     output_nc_file = TRIM(ADJUSTL(outfile_ext))                         &
     &              //'_as_'//sasurf//'_th0_'//stheta0                    &
     &              //'_ph0_'//sphi0//'_ph_'//sphi//'.nc'
  ENDIF

  PRINT *
  PRINT *, '--------------------------------------------------------------------------'
  PRINT *, TRIM(ADJUSTL(output_nc_file))
  PRINT *

  !---------------------------------------------------------------------
  ! Count number of photons that we will use
  !--------------------------------------------------------------------- 
  IF (chatty) PRINT *,'Count number of photons that we will use'

  N_PHOTONS_PER_ROUND = N_PHOTONS_PER_COLUMN * n_columns  !per round
  N_PHOTONS           = N_PHOTONS_PER_ROUND / nt          !per thread

  IF (chatty) PRINT *,'n_columns',n_columns
  IF (chatty) PRINT *,'N_PHOTONS_PER_COLUMN',N_PHOTONS_PER_COLUMN
  IF (chatty) PRINT *,'N_ROUNDS',N_ROUNDS  
  IF (chatty) PRINT *,'N_PHOTONS_PER_ROUND',N_PHOTONS_PER_ROUND
  IF (chatty) PRINT *,'nt',nt
  IF (chatty) PRINT *,'N_PHOTONS',N_PHOTONS

  IF (N_PHOTONS.LT.(nx*ny)) THEN
     CALL write_warning('Less that 1 photon per column used')
     IF (chatty) PRINT *,'N_PHOTONS',N_PHOTONS
     IF (chatty) PRINT *,'nt',nt
     IF (chatty) PRINT *,'nx',nx
     IF (chatty) PRINT *,'ny',ny
     IF (chatty) PRINT *,'n_columns',n_columns
     READ(*,*)
  ENDIF

  !---------------------------------------------------------------------
  ! Swap light source and detector location in case of backward mode
  !---------------------------------------------------------------------  
  IF (chatty) PRINT *,'Swap light source and detector location in case of backward mode'

  IF (backward_mode) THEN
     detector_location_save = detector_location
     detector_location = source_location
     source_location = detector_location_save
  ENDIF

  !---------------------------------------------------------------------
  ! Find k_source (the z layer corresponding to where the light source is)
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Find k_source (the z layer corresponding to where the light source is)'

  IF (source_location .EQ. 1) THEN
     k_source = nz !light source is at TOA
     IF (theta0 .GE. 90.D0) THEN
        PRINT*,'Error: light source is at TOA but sza >= 90 deg'
        STOP
     ENDIF
  ELSE
     k_source = 0 !light source is at bottom layer
     IF (theta0 .LT. 90.D0) THEN
        PRINT *,'theta0=',theta0
        PRINT*,'Error: light source is at bottom layer but sza < 90 deg'
        STOP
     ENDIF
  ENDIF

  !---------------------------------------------------------------------
  ! Find k_detector (the z layer corresponding to where the sensor is)
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Find k_detector (the z layer corresponding to where the sensor is)'

  IF (detector_location .EQ. 1) THEN
     k_detector = nz !detector is at TOA (in vacuum layer)
     signal_name = 'radiance'
     loc_name = 'TOA '
     IF (MAXVAL(vartheta) .GE. 90.D0) THEN
        PRINT*,'Error: detector is at TOA but vza >= 90 deg'
        STOP
     ENDIF
  ELSEIF (detector_location .EQ. 0) THEN
     k_detector = 0  !detector is at surface
     signal_name = 'radiance'
     loc_name = 'Surface '
     IF (MINVAL(vartheta) .LT. 90.D0) THEN
        PRINT*,'Error: detector is at surface but vza < 90 deg'
        STOP
     ENDIF     
  ELSE
     PRINT *,'Error: location of the detector should be TOA (1) or surface (0)'
     STOP
  ENDIF

  !---------------------------------------------------------------------
  ! Find k_register (the z layer where the results are saved)
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Find k_register (the z layer where the results are saved)'

  IF (register_location .EQ. 1) THEN
     k_register = nz !save grid is at TOA (in vacuum layer)
  ELSEIF (register_location .EQ. 0) THEN
     k_register = 0  !save grid is at surface
  ELSE
     PRINT *,'Error: location of the grid where results are saved should be TOA (1) or surface (0)'
     STOP
  ENDIF

  IF (backward_mode) THEN
     IF ((k_source .EQ. 1) .AND. (k_register .NE. 0)) THEN
        PRINT*,'k_register =',k_register
        PRINT*,'k_source=',k_source
        PRINT*,'Error in main.f90: Light source is at surface (k_source .EQ. 1),'
        PRINT*,'but photons are not registered at surface (k_register .NE. 0).'
        PRINT*,'This option is not yet available in MONKI backward mode.'
        STOP
     ENDIF
     IF ((k_source .EQ. nz) .AND. (k_register .NE. nz)) THEN             
        PRINT*,'k_register =',k_register
        PRINT*,'k_source=',k_source
        PRINT*,'Error in main.f90: Light source is at TOA (k_source .EQ. nz),'
        PRINT*,'but photons are not registered at TOA (k_register .NE. nz).'
        PRINT*,'This option is not yet available in MONKI backward mode.'
        STOP
     ENDIF     
  ENDIF

  !---------------------------------------------------------------------
  ! Compute zenith angle of sunlight propagation direction
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Compute zenith angle of sunlight propagation direction'

  vartheta0 = 180.D0 - theta0

  !---------------------------------------------------------------------
  ! Specify surface height and type
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Specify surface height and type'

  i_surface_elev = 0

  DO i=1,nx
     DO j=1,ny
        i_surface_elev(i,j) = 1
        surface(i,j) = -1 ! -1 = Lambertian
     ENDDO
  ENDDO

  !---------------------------------------------------------------------
  ! Populate surface albedo array
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Populate surface albedo array'

  DO i=1,nx
    DO j=1,ny
      IF (surface(i,j)==-1) THEN
        surface_alb_p(i,j) = asurf
      ELSE
        PRINT *,'ERROR: surface is non-Lambertian'
        STOP
      ENDIF
    ENDDO
  ENDDO

  !---------------------------------------------------------------------
  ! Open bmsca and bmabs files if needed
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Open bmsca and bmabs files if needed'

  IF (read_bmsca_from_file) THEN
     bmsca_file = TRIM(ADJUSTL(monki_home))//'/INPUT/'//TRIM(ADJUSTL(bmsca_file))
     OPEN(unit=ibmsca,file=bmsca_file)
     READ(ibmsca,*)
  ENDIF

  IF (read_bmabs_from_file) THEN
     bmabs_file = TRIM(ADJUSTL(monki_home))//'/INPUT/'//TRIM(ADJUSTL(bmabs_file))
     OPEN(unit=ibmabs,file=bmabs_file)
     READ(ibmabs,*)
  ENDIF

  !---------------------------------------------------------------------
  ! Loop over wavelengths or wavenumbers
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Loop over wavelengths or wavenubmers'

  wave = wave_min - dwave
  iw = 0
  nbands = nwaves

10 wave = wave + dwave
  iw = iw + 1
  iband = iw

  !---------------------------------------------------------------------
  ! Get wavelength in nm and wavenumber in cm-1
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Get wavelength in nm and wavenumber in cm-1'

  IF (wave_type .EQ. 'wn') THEN
     wav   = 1.D7/wave !nm
     wavnr = wave      !cm-1
  ELSE 
     wav   = wave      !nm
     wavnr = 1.D7/wave !cm-1
  ENDIF

  PRINT *,'Processing wavelength: ',wav,'nm'

  !---------------------------------------------------------------------
  ! Calculate bmsca and bmabs for each layer
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Calculate bmsca and bmabs for each layer'

  IF (background .EQ. 'Earth') THEN
     CALL molecules_earth(nz,pres,temp,alt,wav,depol,gas_col_type,&
             &           bmsca,monki_home,ngases,vmrgas,bmabs)
  ELSEIF (background .EQ. 'Venus') THEN
     CALL molecules_venus(nz,pres,temp,alt,wav,depol,gas_col_type,&
                         bmsca,monki_home,ngases,vmrgas,bmabs)
  ELSE
     PRINT *,'background = ',background
     PRINT *,'Wrong background specified, should be Earth or Venus.'
     STOP
  ENDIF 

  !---------------------------------------------------------------------
  ! Overwrite depolarization factor if needed
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Overwrite depolarization factor if needed'

  IF (fix_depol) THEN
     depol = depol_fixed
  ENDIF

  !---------------------------------------------------------------------
  ! Overwrite bmsca and/or bmabs if provided in a separate data file
  ! This option is useful for testing and/or code comparisons
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Overwrite bmsca and/or bmabs if provided in a separate data file'

  IF (read_bmsca_from_file) THEN
     wav_read = 0.D0     
     DO WHILE (DABS(wav_read-wav/1000.D0).GT.1.0D-12) 
        READ(ibmsca,*) wav_read, (bmsca(k), k=1,nz-1)
        bmsca(nz) = 0.D0 !vacuum layer
     ENDDO
     IF (DABS(wav_read-wav/1000.D0).GT.1.0D-12) THEN
        PRINT*,'Error: wav not found in bmsca file.'
        STOP
     ENDIF
  ENDIF

  IF (read_bmabs_from_file) THEN
     wav_read = 0.D0     
     DO WHILE (DABS(wav_read-wav/1000.D0).GT.1.0D-12)
        !PRINT*,'wav_read=',wav_read,'; wav/1000=',wav/1000.D0
        READ(ibmabs,*) wav_read, (bmabs(k), k=1,nz-1)
        bmabs(nz) = 0.D0 !vacuum layer
     ENDDO
     IF (DABS(wav_read-wav/1000.D0).GT.1.0D-12) THEN
        PRINT*,'Error: wav not found in bmabs file.'
        STOP
     ENDIF
  ENDIF

  !---------------------------------------------------------------------
  ! Get single scattering matrices and albedo for every aerosol type
  !---------------------------------------------------------------------
  IF (chatty) PRINT*,'Get single scattering matrices and albedo for every aerosol type'

  CALL aerosols(nmat,nscatt_types,new_scatt_info_pol,monki_home,coelut,coefile,&
  &             wav,cloud_input_type,cloud_abundance,nx,ny,nz,baer)

  !---------------------------------------------------------------------
  ! Initialize extinctions and single scattering albedo
  !---------------------------------------------------------------------
  IF (chatty) PRINT*,'Initialize extinctions and single scattering albedo'
  
  BETA=0.D0
  BETA_gas = 0.D0
  ext_by_scatt_species = 0.D0
  norm_ext_by_scatt_species = 0.D0
  ssalb_gas = 1.D0

  !---------------------------------------------------------------------
  ! Define total extinctions (abs+sca) per km and single scat albedo for gas
  !---------------------------------------------------------------------
  IF (chatty) PRINT*,'Define total extinctions (abs+sca) per km and single scat albedo for gas'

  DO k=1,nz

     !PRINT*,'k=',k,'alt(k)=',alt(k),'bmsca(k)=', bmsca(k),' bmabs(k)=', bmabs(k),' z_unit(k)=',z_unit(k)
     BETA_gas(k) = ( bmsca(k) + bmabs(k) ) / z_unit(k) !(bmsca+bmabs)/km

     DO j=1,ny
        DO i=1,nx
           ext_by_scatt_species(nscatt_types+1,i,j,k) =                &
  &                ( bmsca(k) + bmabs(k) ) / z_unit(k) !(bmsca+bmabs)/km

           IF ((bmsca(k) + bmabs(k)) .GT. eps ) THEN
              ssalb_gas(k) = bmsca(k) / ( bmsca(k) + bmabs(k) )
           ENDIF

        ENDDO
     ENDDO
  ENDDO

  !---------------------------------------------------------------------
  ! Define total extinctions (abs+sca) per km for aerosols
  !---------------------------------------------------------------------
  IF (chatty) PRINT*,'Define total extinctions (abs+sca) per km for aerosols'

  DO isc=1,nscatt_types
     DO k=1,nz
        DO j=1,ny
           DO i=1,nx
              ext_by_scatt_species(isc,i,j,k) = baer(isc,i,j,k) / z_unit(k) !(basca+baabs)/km 
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !---------------------------------------------------------------------
  ! Calculate total optical thickness per kilometer (aerosols + gas)
  !---------------------------------------------------------------------
  IF (chatty) PRINT*,'Calculate total optical thickness per kilometer (aerosols + gas)'
  
  DO isc=1,nscatt_types+1
     BETA=BETA+ext_by_scatt_species(isc,:,:,:) !(bsca+babs)/km
  ENDDO

  !---------------------------------------------------------------------
  ! Normalize the component extinctions at each grid point 
  !--------------------------------------------------------------------- 
  IF (chatty) PRINT*,'Normalize the component extinctions at each grid point '

  DO isc=1,nscatt_types+1
     DO i=1,nx
        DO j=1,ny
           DO k=1,nz
              IF (BETA(i,j,k) .GT. eps) THEN
                 norm_ext_by_scatt_species(isc,i,j,k) =                &
 &                    ext_by_scatt_species(isc,i,j,k) / BETA(i,j,k)
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !---------------------------------------------------------------------
  ! Precompute path optical thickness of the gas from bottom of a layer
  ! to the sensor, and per VZA
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Precompute path optical thickness of the gas per z-level and per VZA'

  TAU_gas = 0.D0

  DO iangles=1,n_angles
     DO k0 = 1,nz
        IF (k_detector .EQ. nz) THEN
           DO k=k0,nz
              TAU_gas(k0,iangles) = TAU_gas(k0,iangles) + &
              &   BETA_gas(k)*z_unit(k)/DABS(DCOS(DegtoRad*vartheta(iangles)))
           ENDDO
        ELSEIF (k_detector .EQ. 0) THEN
           DO k=1,k0
              TAU_gas(k0,iangles) = TAU_gas(k0,iangles) + &
              &   BETA_gas(k)*z_unit(k)/DABS(DCOS(DegtoRad*vartheta(iangles)))
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  !---------------------------------------------------------------------
  ! Precompute cells for which the path to the sensor only consists of 
  ! gas.
  ! For now, we take the lowest layer for which the atmosphere is 
  ! horizontally homogeneous up to TOA.
  !---------------------------------------------------------------------  
  IF (chatty) PRINT *,'Precompute cells for which the path to the sensor only consists of gas.'

  aerosol_present = 0
  DO isc=1,nscatt_types
     DO k=1,nz
        DO j=1,ny
           DO i=1,nx 
              IF (baer(isc,i,j,k) .GT. eps) THEN
                 aerosol_present(k) = 1
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  k = nz
  lowest_k_of_only_gas = nz

  IF (chatty) PRINT *,'k=',k
  DO WHILE ((aerosol_present(k) .EQ. 0) .AND. (k .GE. 1)) 
     lowest_k_of_only_gas = k
     !PRINT *,'k=',k
     k = k-1
     IF (k.EQ.0) EXIT
  ENDDO 

  !---------------------------------------------------------------------
  ! Flag grid cells for which the path to the sensor is blocked by a 
  ! cloud or aerosols. 
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Flag grid cells for which the path to the sensor is blocked by a cloud or aerosols.'

  blocked_sensor_path(:,:,:,:) = 0

  precompute_detector_paths = .false.
  IF (precompute_detector_paths) THEN
     CALL flag_blocked_sensor_path(nx,ny,nz,n_angles,nscatt_types,x_unit,y_unit,&
          &             z_unit,vartheta,phi,baer,k_detector,blocked_sensor_path)
  ENDIF

  !---------------------------------------------------------------------
  ! Loop over rounds
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Loop over rounds'

  DO iround = 1,N_ROUNDS

     IF (chatty .OR. basicinfo) PRINT*,' '
     IF (chatty .OR. basicinfo) WRITE(string,'(a,i6,a,i6)') 'Round =', iround,';    max #rounds =',N_rounds
     IF (chatty .OR. basicinfo) CALL Write_info(string)

     !---------------------------------------------------------------------
     ! Initalize arrays
     !---------------------------------------------------------------------      
     IF (chatty) PRINT *,'Initalize arrays'
     
     stokes_p    = 0.D0
     weighted_path_store  = 0.D0
     !weighted_boxpath_store  = 0.D0
     sum_of_weights_store = 0.D0
     CONTRIB_STOKES_3D_p = 0.D0
     SINGSCAT_CONTRIB_STOKES_3D_p = 0.D0
     MULTSCAT_CONTRIB_STOKES_3D_p = 0.D0
     stokes_sing_p = 0.D0
     stokes_mult_p = 0.D0

     !---------------------------------------------------------------------
     ! Run Mont Carlo subroutine divided over threads
     !---------------------------------------------------------------------
     IF (chatty) PRINT *,'Run Mont Carlo subroutine divided over threads'

     !$OMP PARALLEL REDUCTION(+:stokes_p,weighted_path_store,&
     !$OMP& sum_of_weights_store,CONTRIB_STOKES_3D_p,SINGSCAT_CONTRIB_STOKES_3D_p,&
     !$OMP& MULTSCAT_CONTRIB_STOKES_3D_p,Stokes_sing_p,Stokes_mult_p)
           
     CALL monte_carlo(nx, ny, nz, nmat, n_angles, &
                & nx_1, ny_1, nz_1, nmat_1, n_angles_1,&
                & nx_2, ny_2, nz_2, nmat_2, n_angles_2,&
                & nscatt_types, new_scatt_info_pol, &
                & n_photons, &
                & vartheta0,phi0,&
                & vartheta(1:n_angles),phi,& 
                & x_unit, y_unit, z_unit, &
                & surface(1:nx,1:ny),surface_alb_p(1:nx,1:ny),&
                & i_surface_elev(1:nx,1:ny), &
                & BETA,norm_ext_by_scatt_species,ssalb_gas(1:nz),depol, &
                & k_detector,&
                & TAU_gas,lowest_k_of_only_gas,&
                & blocked_sensor_path,precompute_detector_paths,&
                & Stokes_p(1:nmat,1:nx,1:ny,1:n_angles),&
                & save_path,&!nz_box,&
                & weighted_path_store(1:nx_2,1:ny_2,1:nz_2,1:n_angles_2),&
                !& weighted_boxpath_store(1:nx_2,1:ny_2,1:nx_2,1:ny_2,1:nz_box_2,1:n_angles_2),&
                & sum_of_weights_store(1:nx_2,1:ny_2,1:n_angles_2),&
                & save_scats,&
                & CONTRIB_STOKES_3D_p(1:nmat_1,1:nx_1,1:ny_1,0:nz_1,1:n_angles_1),&
                & SINGSCAT_CONTRIB_STOKES_3D_p(1:nmat_1,1:nx_1,1:ny_1,0:nz_1,1:n_angles_1),&
                & MULTSCAT_CONTRIB_STOKES_3D_p(1:nmat_1,1:nx_1,1:ny_1,0:nz_1,1:n_angles_1),&
                & Stokes_sing_p(1:nmat_1,1:nx_1,1:ny_1,1:n_angles_1),&
                & Stokes_mult_p(1:nmat_1,1:nx_1,1:ny_1,1:n_angles_1),&
                & k_register,k_source,backward_mode)

     !$OMP END PARALLEL

     !---------------------------------------------------------------------
     ! Accumulate results for this round
     !---------------------------------------------------------------------
     IF (chatty) PRINT *,'Accumulate results for this round'

     CALL accumulate(iround,finished_rad,tol_error)

     !---------------------------------------------------------------------
     ! Exit loop in rounds when accuracy is sufficient
     !---------------------------------------------------------------------
     IF (finished_rad .AND. stop_converged) THEN
        IF (chatty .OR. basicinfo) PRINT*,'Result converged.'
        EXIT
     ENDIF
          
  !---------------------------------------------------------------------
  ! End loop over rounds
  !---------------------------------------------------------------------
  ENDDO 

  IF (chatty) PRINT *,'End loop over rounds'

  !---------------------------------------------------------------------
  ! After EXIT DO loop, FORTRAN increments i_round by 1! 
  ! Set i_round back to N_rounds:     
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Set i_round back to N_rounds:'

  IF (iround .GT. N_rounds) THEN
     iround = N_rounds
  ENDIF

  !---------------------------------------------------------------------
  ! Divide fluxes and radiances by #rounds (N) and N_PHOTONS_PER_COLUMN
  !---------------------------------------------------------------------                      
  IF (chatty) PRINT *,'Divide fluxes and radiances by #rounds (N) and N_PHOTONS_PER_COLUMN'

  SUM_STOKES = SUM_STOKES/iround/N_PHOTONS_PER_COLUMN      
  
  IF (save_scats) THEN

     SUM_STOKES_SING = SUM_STOKES_SING/iround/N_PHOTONS_PER_COLUMN
     SUM_STOKES_MULT = SUM_STOKES_MULT/iround/N_PHOTONS_PER_COLUMN

     SUM_CONTRIB_STOKES_3D = SUM_CONTRIB_STOKES_3D/iround/N_PHOTONS_PER_COLUMN
     SUM_SINGSCAT_CONTRIB_STOKES_3D = SUM_SINGSCAT_CONTRIB_STOKES_3D/iround/N_PHOTONS_PER_COLUMN
     SUM_MULTSCAT_CONTRIB_STOKES_3D = SUM_MULTSCAT_CONTRIB_STOKES_3D/iround/N_PHOTONS_PER_COLUMN
  
  ENDIF

  !---------------------------------------------------------------------
  ! Reverse sign of U in case of backward Monte Carlo
  !---------------------------------------------------------------------   
  IF (backward_mode) THEN
     IF (nmat .GE. 3) THEN             

        IF (chatty) PRINT *,'Reverse sign of U in case of backward Monte Carlo'     
        
        SUM_STOKES(3,:,:,:) = -SUM_STOKES(3,:,:,:)
        
        IF (save_convergence_history) THEN
           STOKES_ROUNDS_save(3,:) = -STOKES_ROUNDS_save(3,:)
        ENDIF

        IF (save_scats) THEN
        
           SUM_STOKES_SING(3,:,:,:) = -SUM_STOKES_SING(3,:,:,:)
           SUM_STOKES_MULT(3,:,:,:) = -SUM_STOKES_MULT(3,:,:,:)
        
           SUM_CONTRIB_STOKES_3D(3,:,:,:,:) = -SUM_CONTRIB_STOKES_3D(3,:,:,:,:)
           SUM_SINGSCAT_CONTRIB_STOKES_3D(3,:,:,:,:) = -SUM_SINGSCAT_CONTRIB_STOKES_3D(3,:,:,:,:)
           SUM_MULTSCAT_CONTRIB_STOKES_3D(3,:,:,:,:) = -SUM_MULTSCAT_CONTRIB_STOKES_3D(3,:,:,:,:) 

        ENDIF
     ENDIF
  ENDIF

  !---------------------------------------------------------------------
  ! Compute degree of polarization (DOP)
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Compute degree of polarization (DOP)'

  DOP = 0.D0

  DO i=1,nx
     DO j=1,ny
        DO ia=1,n_angles
           IF (SUM_STOKES(1,i,j,ia) .GT. eps) THEN
              IF (nmat .EQ. 3) THEN
                 DOP(i,j,ia) = DSQRT( SUM_STOKES(2,i,j,ia)**2   &
                           & + SUM_STOKES(3,i,j,ia)**2 ) &
                           & / SUM_STOKES(1,i,j,ia)
              ELSE IF (nmat .EQ. 4) THEN 
                 DOP(i,j,ia) = DSQRT( SUM_STOKES(2,i,j,ia)**2   &
                           & + SUM_STOKES(3,i,j,ia)**2   &
                           & + SUM_STOKES(4,i,j,ia)**2 ) &
                           & / SUM_STOKES(1,i,j,ia)
              ENDIF
           ENDIF   
        ENDDO
     ENDDO
  ENDDO

  !---------------------------------------------------------------------
  ! For each pixel (i,j,ia), compute the layer (k) dependent mean light
  ! paths and air mass factors (AMFs)
  !---------------------------------------------------------------------
  IF (save_path) THEN

     IF (chatty) PRINT*,'For each pixel (i,j,ia), compute the layer (k) dependent mean light path and AMFs'
  
     DO k=1,nz
        mean_path_layer(:,:,k,:) = SUM_weighted_path_store(:,:,k,:) / SUM_sum_of_weights_store(:,:,:)
        AMF_layer(:,:,k,:) = mean_path_layer(:,:,k,:) / z_unit(k)

     ENDDO

  ENDIF

  !IF (save_path) THEN
  !              
  !   DO k=1,nz_box
  !      DO j=1,ny
  !         DO i=1,nx  
  !            mean_path_box(:,:,i,j,k,:) = SUM_weighted_boxpath_store(:,:,i,j,k,:) / SUM_sum_of_weights_store(:,:,:)
  !            AMF_box(:,:,i,j,k,:) = mean_path_box(:,:,i,j,k,:) / z_unit(k)
  !         ENDDO
  !      ENDDO
  !   ENDDO  
  !
  !ENDIF

  !---------------------------------------------------------------------
  ! Write results for this band to output file
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Writing data to file for current band'

  CALL Report_ncdf(iround,iband,output_nc_file)

  IF (chatty) PRINT *,'Finished Writing data to file'

  !---------------------------------------------------------------------
  ! Next wavelength or end loop over wavelengths 
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Next wavelength or end loop over wavelengths'
  
  IF ((wave + dwave) .LE. (wave_max)) GOTO 10

  !---------------------------------------------------------------------
  ! Close bmsca and bmabs files if needed
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Close bmsca and bmabs files if needed'

  IF (read_bmsca_from_file) THEN
     CLOSE(unit=ibmsca)
  ENDIF

  IF (read_bmabs_from_file) THEN
     CLOSE(unit=ibmabs)
  ENDIF

  !---------------------------------------------------------------------
  ! Get the end time
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Get the end time'
  
  CALL system_clock(end_clock, clock_rate)
  IF (chatty) PRINT*,'End clock: ',end_clock

  !---------------------------------------------------------------------
  ! Calculate and print run time
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Calculate and print run time'

  ! Calculate elapsed time in seconds
  elapsed_seconds = REAL(end_clock - start_clock) / REAL(clock_rate)

  ! Convert elapsed time to total seconds
  total_seconds = INT(elapsed_seconds)
  milliseconds = (elapsed_seconds - total_seconds) * 1000.0

  ! Calculate days, hours, minutes, and seconds
  days    = total_seconds / (24 * 3600)
  hours   = MOD(total_seconds, 24 * 3600) / 3600
  minutes = MOD(total_seconds, 3600) / 60
  seconds = MOD(total_seconds, 60)

  ! Print the elapsed time
  print *, "Elapsed time: ", days, " days, ", hours, " hours, ", minutes, " minutes, ", &
             seconds, " seconds, ", milliseconds, " milliseconds."

  !---------------------------------------------------------------------
  ! Finish
  !---------------------------------------------------------------------
  IF (chatty) PRINT *,'Finish'

  IF (status .NE. 0) THEN
     CALL write_error(error_str)
     CALL write_error('An error has occured.')
     CALL EXIT(1)
  ELSE
     CALL write_info('MONKI FINISHED with no reported errors.')
     CALL EXIT(0)
  ENDIF

  !---------------------------------------------------------------------
  !**************************** END SCRIPT *****************************
  !---------------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------------
!********************** INTERNAL SUBROUTINES *************************
!---------------------------------------------------------------------

! NOTE THAT VARIABLES ARE SHARED WITH PROGRAM ABOVE!
! That is, passing input and output arguments explicitly is not 
! required, but it is good practice for code readability.

!+-----------------------------------------------------------------
! SUBROUTINE accumulate
!+----------------------------------------------------------------- 
! Accummulates results after a each round, and computes the 
! error. If the error is small enough, finshed_rad is set equal
! to .true.
!+-----------------------------------------------------------------
  SUBROUTINE accumulate(iround,finished_rad,tol_error) 
    
    ! Input --------------------------------------------------------------
    INTEGER,INTENT(in)    :: iround    !Current round number
    REAL*8,INTENT(in)     :: tol_error !Tolerance error for convergence

    ! In- and output -----------------------------------------------------
    LOGICAL,INTENT(inout) :: finished_rad !Indicates if we are finished

    !---------------------------------------------------------------------
    ! Collect results
    !---------------------------------------------------------------------
    STOKES = STOKES_p(1:nmat,1:nx,1:ny,1:n_angles)
    
    IF (save_scats) THEN
       STOKES_SING = STOKES_SING_p(1:nmat,1:nx,1:ny,1:n_angles)
       STOKES_MULT = STOKES_MULT_p(1:nmat,1:nx,1:ny,1:n_angles)
       CONTRIB_STOKES_3D = CONTRIB_STOKES_3D_p(1:nmat,1:nx,1:ny,0:nz,1:n_angles)
       SINGSCAT_CONTRIB_STOKES_3D = SINGSCAT_CONTRIB_STOKES_3D_p(1:nmat,1:nx,1:ny,0:nz,1:n_angles)
       MULTSCAT_CONTRIB_STOKES_3D = MULTSCAT_CONTRIB_STOKES_3D_p(1:nmat,1:nx,1:ny,0:nz,1:n_angles)
    ENDIF

    !---------------------------------------------------------------------
    ! In case of first round (iround = current round nr):
    !---------------------------------------------------------------------
    IF (iround == 1) THEN
       STD_STOKES = 999.D0
       SUM_STOKES = 0.D0
       SUM_STOKES_ROUNDS = 0.D0
       SUM_STOKES_ROUNDS_SQ = 0.D0
       finished_rad = .False.
       ERROR_STOKES = 0.D0
       IF (save_convergence_history) THEN
          !S1_ROUNDS_save = 0.D0
          STOKES_ROUNDS_save = 0.D0
       ENDIF
       IF (save_path) THEN
          SUM_weighted_path_store = 0.D0
          !SUM_weighted_boxpath_store = 0.D0
          SUM_sum_of_weights_store = 0.D0
       ENDIF
       IF (save_scats) THEN
          SUM_STOKES_SING = 0.D0
          SUM_STOKES_MULT = 0.D0 
          SUM_CONTRIB_STOKES_3D = 0.D0    
          SUM_SINGSCAT_CONTRIB_STOKES_3D = 0.D0
          SUM_MULTSCAT_CONTRIB_STOKES_3D = 0.D0
       ENDIF
    ENDIF

    !---------------------------------------------------------------------
    ! Update the sums and errors
    !---------------------------------------------------------------------       
    SUM_STOKES    = SUM_STOKES + STOKES

    IF (save_path) THEN
       SUM_sum_of_weights_store = SUM_sum_of_weights_store + sum_of_weights_store(1:nx,1:ny,1:n_angles)
    ENDIF

    IF (save_path) THEN
       SUM_weighted_path_store = SUM_weighted_path_store + weighted_path_store(1:nx,1:ny,1:nz,1:n_angles)
    ENDIF

    !IF (save_path) THEN
    !   SUM_weighted_boxpath_store = SUM_weighted_boxpath_store + weighted_boxpath_store(1:nx,1:ny,1:nx,1:ny,1:nz_box,1:n_angles)
    !ENDIF

    IF (save_scats) THEN
       SUM_STOKES_SING = SUM_STOKES_SING + STOKES_SING
       SUM_STOKES_MULT = SUM_STOKES_MULT + STOKES_MULT

       SUM_CONTRIB_STOKES_3D = SUM_CONTRIB_STOKES_3D + &
               & CONTRIB_STOKES_3D(1:nmat,1:nx,1:ny,0:nz,1:n_angles)
       SUM_SINGSCAT_CONTRIB_STOKES_3D = SUM_SINGSCAT_CONTRIB_STOKES_3D + &
               & SINGSCAT_CONTRIB_STOKES_3D(1:nmat,1:nx,1:ny,0:nz,1:n_angles)
       SUM_MULTSCAT_CONTRIB_STOKES_3D = SUM_MULTSCAT_CONTRIB_STOKES_3D + &
               & MULTSCAT_CONTRIB_STOKES_3D(1:nmat,1:nx,1:ny,0:nz,1:n_angles)
    ENDIF

    !---------------------------------------------------------------------
    ! Accumulate results of different rounds
    !---------------------------------------------------------------------
    SUM_STOKES_ROUNDS   = SUM_STOKES_ROUNDS + STOKES/N_PHOTONS_PER_COLUMN
    SUM_STOKES_ROUNDS_SQ = SUM_STOKES_ROUNDS_SQ + (STOKES/N_PHOTONS_PER_COLUMN)**2

    !---------------------------------------------------------------------
    ! Save estimate per round
    !---------------------------------------------------------------------
    IF (save_convergence_history) THEN
       !S1_ROUNDS_save(iround) = SUM(SUM_STOKES_ROUNDS(1,:,:,:))/nx/ny/n_angles/REAL(iround,8)
       STOKES_ROUNDS_save(:,iround) = SUM(SUM(SUM(SUM_STOKES_ROUNDS(:,:,:,:),dim=2),dim=2),dim=2)/nx/ny/n_angles/REAL(iround,8)
    ENDIF

    !---------------------------------------------------------------------
    ! Compute standard deviation for the 30th round and higher
    !---------------------------------------------------------------------    
    IF (iround .GT. 30) THEN

       STD_STOKES = DSQRT( (SUM_STOKES_ROUNDS_SQ - 1.D0/iround * SUM_STOKES_ROUNDS**2) / (iround-1) ) / DSQRT(REAL(iround,8))

       !---------------------------------------------------------------------
       ! IF the max error of the radiance in the column < tol_error THEN we
       ! are done.
       !---------------------------------------------------------------------
       STD_S1_MAX = MAXVAL(STD_STOKES(1,:,:,:))
       STD_S1_MEAN = SUM(STD_STOKES(1,:,:,:))/nx/ny/n_angles
       
       IF ( (STD_S1_MAX.LE.tol_error) .AND. (iround .GE. 30) ) THEN        
          finished_rad=.True.
       ELSE
          finished_rad=.False.
       ENDIF

       !---------------------------------------------------------------------
       ! Print the results to screen.
       !---------------------------------------------------------------------      
       IF (chatty .OR. basicinfo) PRINT*,'Mean '//TRIM(ADJUSTL(signal_name))//' std  = ',STD_S1_MEAN,&
               &       ', Max '//TRIM(ADJUSTL(signal_name))//' std  = ',STD_S1_MAX,', Goal = ',tol_error

 
    ENDIF

       
    RETURN

  END SUBROUTINE accumulate


!+-----------------------------------------------------------------
! SUBROUTINE report_ncdf
!+-----------------------------------------------------------------
! Writes output to the ncdf output file
!+-----------------------------------------------------------------
  SUBROUTINE report_ncdf(iround,iband,outfilename)
    
    ! Input -------------------------------------------------------      
    INTEGER,INTENT(in) :: iround !Current round number
    INTEGER,INTENT(in) :: iband  !Current wavelength number
    CHARACTER(len=500),INTENT(in) :: outfilename !output filename

    ! Internal ----------------------------------------------------
    INTEGER,SAVE  :: xdim,ydim,zdim,levdim,gasdim,ang_dim,matDim!zboxdim,n_sites
    INTEGER,SAVE  :: ncid,bands_dim,single_dim,double_dim,zdim_withsurf
    INTEGER,SAVE  :: scattdim,nroundsDIM
    INTEGER,SAVE  :: nz_novac,nlev_novac,nz_withsurf_novac

    INTEGER,SAVE :: x_unitID,y_unitID,z_unitID
    INTEGER,SAVE :: vmrgasID,presID,tempID,altID
    INTEGER,SAVE :: asurfID,reflconvID
    INTEGER,SAVE :: theta0ID,phi0ID,phiID,thetaID
    INTEGER,SAVE :: N_PHOTONS_PER_COLUMN_ID,N_ROUNDS_ID


    CHARACTER(len=500)    :: err_str,outfile
    
    !---------------------------------------------------------------------
    ! If outputnum == 0 then we create a new ouput file, define all the 
    ! variables and write out the general information
    !---------------------------------------------------------------------
    IF (outputnum.EQ.0) THEN
    
       !---------------------------------------------------------------------
       ! Remove vacuum layer
       !---------------------------------------------------------------------
       nz_novac = nz-1
       nlev_novac = nlev-1
       nz_withsurf_novac = nz

       !---------------------------------------------------------------------
       ! Create output file, in mode nf90_netcdf4
       ! ncid is the ID for this file in what follows
       !---------------------------------------------------------------------
       outfile = TRIM(ADJUSTL(monki_home))//'/OUTPUT/'//TRIM(ADJUSTL(outfilename))
       IF (chatty) write(*,*) 'Initalizing:',outfile
       err_str=' Error in nf90_create'
       status = nf90_create(outfile, nf90_netcdf4, ncid)

       !---------------------------------------------------------------------
       ! Define dimensions
       !---------------------------------------------------------------------
       IF (chatty) PRINT *,'Define dimensions'
       CALL my_nf90_def_dim(ncid, "nx", nx, xDim, status)
       err_str='def_dim_nx'
       IF (status /= 0) GOTO 999
       CALL my_nf90_def_dim(ncid, "ny", ny, yDim, status)
       err_str='def_dim_ny'
       IF (status /= 0) GOTO 999
       CALL my_nf90_def_dim(ncid, "nz", nz_novac, zDim, status)
       err_str='def_dim_nz'
       CALL my_nf90_def_dim(ncid, "nz_withsurf", nz_withsurf_novac, zDim_withsurf, status)
       err_str='def_dim_nz_withsurf'       
       !CALL my_nf90_def_dim(ncid, "nz_box", nz_box, zboxDim, status)
       !err_str='def_dim_nz_box'
       CALL my_nf90_def_dim(ncid, "nlev", nlev_novac, levDim, status)
       err_str='def_dim_nlev'       
       CALL my_nf90_def_dim(ncid, "ngases", ngases, gasDim, status)
       err_str='def_dim_ngases' 
       CALL my_nf90_def_dim(ncid, "nscatt_types", nscatt_types, scattDim, status)
       err_str='def_dim_nscatt_types'      
       CALL my_nf90_def_dim(ncid, "nmat", nmat, matDim, status)
       err_str='def_dim_nmat'       
       IF (status /= 0) GOTO 999
       IF (backward_mode) THEN
          CALL my_nf90_def_dim(ncid, "nsza", n_angles, ang_dim, status)
          err_str='def_dim_nangles'
       ELSE
          CALL my_nf90_def_dim(ncid, "nvza", n_angles, ang_dim, status)
          err_str='def_dim_nangles'          
       ENDIF
       IF (status /= 0) GOTO 999
       CALL my_nf90_def_dim(ncid, "single",  1, single_dim,status)
       err_str='def_dim_single'
       IF (status /= 0) GOTO 999
       CALL my_nf90_def_dim(ncid, "double",  2, double_dim,status)
       err_str='def_dim_double'
       IF (status /= 0) GOTO 999
       CALL my_nf90_def_dim(ncid, "nwavs",  nbands, bands_dim, status)
       err_str='def_dim_nbands'
       IF (status /= 0) GOTO 999
       CALL my_nf90_def_dim(ncid, "nrounds", n_rounds, nroundsDim, status)
       err_str='def_dim_nrounds'
       IF (status /= 0) GOTO 999       

       IF (chatty) PRINT *,'Define variables'

       !----------------------
       name='asurf'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/single_dim/), asurfID, status)
       If (status /= 0) Goto 999         

       !----------------------
       !name='nx'
       !err_str='def_var'//name
       !CALL my_nf90_def_var(ncid, name,  nf90_int, (/single_dim/), nxID, status)
       !If (status /= 0) Goto 999

       !name='ny'
       !err_str='def_var'//name
       !CALL my_nf90_def_var(ncid, name,  nf90_int, (/single_dim/), nyID, status)
       !If (status /= 0) Goto 999

       !name='nz'
       !err_str='def_var'//name
       !CALL my_nf90_def_var(ncid, name,  nf90_int, (/single_dim/), nzID, status)
       !If (status /= 0) Goto 999

       name='n_photons_per_column'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_int, (/single_dim/), n_photons_per_column_ID, status)
       If (status /= 0) Goto 999

       name='n_rounds'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_int, (/single_dim/), n_rounds_ID, status)
       If (status /= 0) Goto 999         

       ! ---------------------

       IF (backward_mode) THEN
          name='theta0'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/ang_dim/), theta0ID, status)
          If (status /= 0) Goto 999

          name='theta'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/single_dim/), thetaID, status)
          If (status /= 0) Goto 999
       ELSE
          name='theta0'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/single_dim/), theta0ID, status)
          If (status /= 0) Goto 999

          name='theta'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/ang_dim/), thetaID, status)
          If (status /= 0) Goto 999
       ENDIF

       name='phi0'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/single_dim/), phi0ID, status)
       If (status /= 0) Goto 999

       name='phi'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/single_dim/), phiID, status)
       If (status /= 0) Goto 999


       ! ---------------------
       name='x_unit'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/single_dim/), x_unitID, status)
       If (status /= 0) Goto 999         

       name='y_unit'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/single_dim/), y_unitID, status)
       If (status /= 0) Goto 999

       name='z_unit'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/zDim/), z_unitID, status)
       If (status /= 0) Goto 999         

       ! ---------------------
       name='vmrgas'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/gasDim,levDim/), vmrgasID, status)
       If (status /= 0) Goto 999

       name='pres'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/levDim/), presID, status)
       If (status /= 0) Goto 999         
      
       name='temp'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/levDim/), tempID, status)
       If (status /= 0) Goto 999

       name='alt'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/levDim/), altID, status)
       If (status /= 0) Goto 999

       name='bmabs'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,zDim/), bmabsID, status)
       If (status /= 0) Goto 999

       name='bmsca'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,zDim/), bmscaID, status)
       If (status /= 0) Goto 999

       name='baer'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,scattDim,xDim,yDim,zDim/), baerID, status)
       If (status /= 0) Goto 999

       name='wav'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim/), wavID, status)


       ! ---------------------
       signal_name = TRIM(ADJUSTL(signal_name))
       loc_with_signal_name =  TRIM(ADJUSTL(loc_name))//' '//TRIM(ADJUSTL(signal_name))

       name=signal_name
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,matDim, xDim, yDim, ang_dim/), stokesID, status)
       If (status /= 0) Goto 999

       IF (save_convergence_history) THEN
          name='radiance_convergence'
          err_str='def_var'//name
          !CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,nroundsDim/), reflconvID, status)
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,matDim,nroundsDim/), reflconvID, status)
          If (status /= 0) Goto 999
       ENDIF

       name='std_'//signal_name
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,matDim, xDim, yDim, ang_dim/), stdstokesID, status)
       If (status /= 0) Goto 999       

       IF (save_scats) THEN

          name= TRIM(ADJUSTL(signal_name))//'_sing'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,matDim, xDim, yDim, ang_dim/), stokessingID, status)
          If (status /= 0) Goto 999

          name= TRIM(ADJUSTL(signal_name))//'_mult'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,matDim, xDim, yDim, ang_dim/), stokesmultID, status)
          If (status /= 0) Goto 999

          name='contrib_'// TRIM(ADJUSTL(signal_name))//'_3D'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim, matDim, xDim, yDim, zDim_withsurf, ang_dim/), &
                  & contrib_stokesID, status)

          name='contrib_'// TRIM(ADJUSTL(signal_name))//'_3D_sing'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,matDim, xDim, yDim, zDim_withsurf, ang_dim/), &
                    & singscat_contrib_stokesID, status)
          If (status /= 0) Goto 999         

          name='contrib_'// TRIM(ADJUSTL(signal_name))//'_3D_mult'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,matDim, xDim, yDim, zDim_withsurf, ang_dim/), &
                    & multscat_contrib_stokesID, status)
          If (status /= 0) Goto 999         
       ENDIF

       name='DOP'
       err_str='def_var'//name
       CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,xDim, yDim, ang_dim/), dopID, status)
       If (status /= 0) Goto 999         

       IF (save_path) THEN
          name='AMF_layer'
          err_str='def_var'//name
          CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,xDim, yDim, zDim, ang_dim/), amfID, status)
          If (status /= 0) Goto 999
         
          !name='AMF_box'
          !err_str='def_var'//name
          !CALL my_nf90_def_var(ncid, name,  nf90_double, (/bands_dim,xDim, yDim, xDim, yDim, zboxDim, ang_dim/), boxamfID, status)
          !If (status /= 0) Goto 999
       ENDIF

       !---------------------------------------------------------------------
       ! Add attributes
       !---------------------------------------------------------------------
       CALL add_attributes(ncid,bmabsID,'Molecular absorption optical thickness','-','bmabs' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,bmscaID,'Molecular scattering optical thickness','-','bmsca' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,baerID,'Aerosol total optical thickness (baerabs+baersca)','-','baer' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,wavID,'Wavelength of the light','nm','wav' ,status)             
       If (status /= 0) Goto 999

       CALL add_attributes(ncid,asurfID,'Surface albedo','-','asurf' ,status)
       If (status /= 0) Goto 999

       !CALL add_attributes(ncid,nxID,'Number of cells in the x-direction','-','nx' ,status)
       !If (status /= 0) Goto 999
       !CALL add_attributes(ncid,nyID,'Number of cells in the y-direction','-','ny' ,status)
       !If (status /= 0) Goto 999
       !CALL add_attributes(ncid,nzID,'Number of cells in the z-direction','-','nz' ,status)
       !If (status /= 0) Goto 999
       CALL add_attributes(ncid,n_photons_per_column_ID,&
               &'Number of photons used per column and per round','-','n_photons_per_column' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,n_rounds_ID,'Number of rounds used','-','n_rounds' ,status)
       If (status /= 0) Goto 999

       CALL add_attributes(ncid,x_unitID,'Spatial resolution in the x-direction','km','x_unit' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,y_unitID,'Spatial resolution in the y-direction','km','y_unit' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,z_unitID,'Spatial resolution in the z-direction (layer dependent)','km','z_unit' ,status)
       If (status /= 0) Goto 999

       CALL add_attributes(ncid,vmrgasID,'Volume mixing ratio of the absorbing gases','ppm','vmrgas' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,presID,'Pressure','Pa','pres' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,tempID,'Temperature','K','temp' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,altID,'Altitude','m','alt' ,status)
       If (status /= 0) Goto 999

       CALL add_attributes(ncid,theta0ID,'Solar zenith angle','deg','solar_zenith_angle' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,phi0ID,'Solar azimuth angle','deg','solar_azimuth_angle' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,phiID,'Viewing azimuth angle','deg','viewing_azimuth_angle' ,status)
       If (status /= 0) Goto 999
       CALL add_attributes(ncid,thetaID,'Viewing zenith angle','deg','viewing_zenith_angle' ,status)
       If (status /= 0) Goto 999

       CALL add_attributes(ncid,stokesID,TRIM(ADJUSTL(loc_with_signal_name)),'-, assuming E0=pi',signal_name,status)
       If (status /= 0) Goto 999

       IF (save_convergence_history) THEN
          CALL add_attributes(ncid,reflconvID,'radiance_convergence','-','radiance convergence history' ,status)
          If (status /= 0) Goto 999       
       ENDIF

       CALL add_attributes(ncid,stdstokesID,'Standard deviation of '//TRIM(ADJUSTL(loc_with_signal_name)),'-','std' ,status)
       If (status /= 0) Goto 999

       IF (save_scats) THEN
          CALL add_attributes(ncid,stokessingID,TRIM(ADJUSTL(loc_with_signal_name))//' single scattering',&
                  &'-', TRIM(ADJUSTL(signal_name))//'_sing' ,status)
          If (status /= 0) Goto 999

          CALL add_attributes(ncid,stokesmultID,TRIM(ADJUSTL(loc_with_signal_name))//' multiple scattering',&
                  &'-', TRIM(ADJUSTL(signal_name))//'_mult' ,status)
          If (status /= 0) Goto 999                

          CALL add_attributes(ncid,contrib_stokesID,TRIM(ADJUSTL(loc_with_signal_name))//' 3D contribution',&
                  & '[-]','contrib_'// TRIM(ADJUSTL(signal_name)),status)
          If (status /= 0) Goto 999
          CALL add_attributes(ncid,singscat_contrib_stokesID,TRIM(ADJUSTL(loc_with_signal_name))//&
                  &' 3D singscat contribution','-','contrib_'// TRIM(ADJUSTL(signal_name))//'_sing' ,status)
          If (status /= 0) Goto 999             
          CALL add_attributes(ncid,multscat_contrib_stokesID,TRIM(ADJUSTL(loc_with_signal_name))&
                  &//' 3D multscat contribution','-','contrib_'// TRIM(ADJUSTL(signal_name))//'_mult' ,status)
          If (status /= 0) Goto 999
       ENDIF

       CALL add_attributes(ncid,dopID,'Degree of polarization','-','DOP' ,status)
       If (status /= 0) Goto 999             
     
       IF (save_path) THEN
          CALL add_attributes(ncid,amfID,'Layer AMF','-','AMF_layer' ,status)
          If (status /= 0) Goto 999

          !CALL add_attributes(ncid,boxamfID,'Box AMF','[-]','AMF_box' ,status)
          !If (status /= 0) Goto 999
       ENDIF             
    
    ENDIF

    !---------------------------------------------------------------------
    ! Write variable content to output file
    !---------------------------------------------------------------------
    IF (chatty) PRINT *,'Write variable content to output file'

    status = nf90_put_var(ncid,bmabsID,bmabs(1:nz_novac),start=(/iband,1/),count=(/1,nz_novac/))
    IF (status /= 0) THEN
       CALL write_error('Error in nf90_put_var bmabs')
       GOTO 999
    ENDIF

    status = nf90_put_var(ncid,bmscaID,bmsca(1:nz_novac),start=(/iband,1/),count=(/1,nz_novac/))
    IF (status /= 0) THEN
       CALL write_error('Error in nf90_put_var bmsca')
       GOTO 999
    ENDIF

    status = nf90_put_var(ncid,baerID,baer(:,:,:,1:nz_novac),start=(/iband,1,1,1,1/),count=(/1,nscatt_types,nx,ny,nz_novac/))
    IF (status /= 0) THEN
       CALL write_error('Error in nf90_put_var baer')
       GOTO 999
    ENDIF

    status = nf90_put_var(ncid,wavID,wav,start=(/iband/))
    IF (status /= 0) THEN
       CALL write_error('Error in nf90_put_var wav')
       GOTO 999
    ENDIF

    IF (outputnum.EQ.0) THEN
    
       outputnum = 1
       status = nf90_put_var(ncid,asurfID,asurf)
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var asurf')
          GOTO 999
       ENDIF

       ! ---------------------
       !status = nf90_put_var(ncid,nxID,nx)
       !IF (status /= 0) THEN
       !   CALL write_error('Error in nf90_put_brdf')
       !   GOTO 999
       !ENDIF

       !status = nf90_put_var(ncid,nyID,ny)
       !IF (status /= 0) THEN
       !   CALL write_error('Error in nf90_put_brdf')
       !   GOTO 999
       !ENDIF

       !status = nf90_put_var(ncid,nzID,nz)
       !IF (status /= 0) THEN
       !   CALL write_error('Error in nf90_put_brdf')
       !   GOTO 999
       !ENDIF

       status = nf90_put_var(ncid,n_photons_per_column_ID,n_photons_per_column)
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var n_photons_per_column')
          GOTO 999
       ENDIF

       status = nf90_put_var(ncid,n_rounds_ID,n_rounds)
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var n_rounds')
          GOTO 999
       ENDIF

       ! ---------------------
       IF (backward_mode) THEN

          status = nf90_put_var(ncid,thetaID,theta0)
          IF (status /= 0) THEN
             CALL write_error('Error in nf90_put_var theta0')
             GOTO 999
          ENDIF

          status = nf90_put_var(ncid,theta0ID,vartheta(:),start=(/1/),count=(/n_angles/))
          IF (status /= 0) THEN
             CALL write_error('Error in nf90_put_var vartheta')
             GOTO 999
          ENDIF

          status = nf90_put_var(ncid,phiID,MOD(phi0+180.D0,360.D0))
          IF (status /= 0) THEN
             CALL write_error('Error in nf90_put_var phi0')
             GOTO 999
          ENDIF

          status = nf90_put_var(ncid,phi0ID,MOD(phi+180.D0,360.D0))
          IF (status /= 0) THEN
             CALL write_error('Error in nf90_put_var phi')
             GOTO 999
          ENDIF          

       ELSE

          status = nf90_put_var(ncid,theta0ID,theta0)
          IF (status /= 0) THEN
             CALL write_error('Error in nf90_put_var theta0')
             GOTO 999
          ENDIF

          status = nf90_put_var(ncid,thetaID,vartheta(:),start=(/1/),count=(/n_angles/))
          IF (status /= 0) THEN
             CALL write_error('Error in nf90_put_var vartheta')
             GOTO 999
          ENDIF

          status = nf90_put_var(ncid,phi0ID,phi0)
          IF (status /= 0) THEN
             CALL write_error('Error in nf90_put_var phi0')
             GOTO 999
          ENDIF

          status = nf90_put_var(ncid,phiID,phi)
          IF (status /= 0) THEN
             CALL write_error('Error in nf90_put_var phi')
             GOTO 999
          ENDIF

       ENDIF

       ! ---------------------
       status = nf90_put_var(ncid,x_unitID,x_unit)
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var x_unit')
          GOTO 999
       ENDIF

       status = nf90_put_var(ncid,y_unitID,y_unit)
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var y_unit')
          GOTO 999
       ENDIF       

       status = nf90_put_var(ncid,z_unitID,z_unit(1:nz_novac),start=(/1/),count=(/nz_novac/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var z_unit')
          GOTO 999
       ENDIF       

       ! ---------------------
       status = nf90_put_var(ncid,vmrgasID,vmrgas(:,1:nlev_novac),start=(/1,1/),count=(/ngases,nlev_novac/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var vmrgas')
          GOTO 999
       ENDIF

       status = nf90_put_var(ncid,presID,pres(1:nlev_novac),start=(/1/),count=(/nlev_novac/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var pres')
          GOTO 999
       ENDIF

       status = nf90_put_var(ncid,tempID,temp(1:nlev_novac),start=(/1/),count=(/nlev_novac/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var temp')
          GOTO 999
       ENDIF

       status = nf90_put_var(ncid,altID,alt(1:nlev_novac),start=(/1/),count=(/nlev_novac/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var alt')
          GOTO 999
       ENDIF


    ENDIF

    ! ---------------------
    status = nf90_put_var(ncid,stokesID,SUM_STOKES(:,:,:,:),start=(/ iband,1,1,1,1 /),count=(/1,nmat,nx,ny,n_angles/))
    IF (status /= 0) THEN
       CALL write_error('Error in nf90_put_var SUM_STOKES')
       GOTO 999
    ENDIF       

    IF (save_convergence_history) THEN
       !status = nf90_put_var(ncid,reflconvID,S1_ROUNDS_save(:),start=(/iband,1/),count=(/1,n_rounds/))
       status = nf90_put_var(ncid,reflconvID,STOKES_ROUNDS_save(:,:),start=(/iband,1,1/),count=(/1,nmat,n_rounds/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var STOKES_ROUNDS_save')
          GOTO 999
       ENDIF
    ENDIF

    status = nf90_put_var(ncid,stdstokesID,STD_STOKES(:,:,:,:),start=(/ iband,1,1,1,1 /),count=(/1,nmat,nx,ny,n_angles/))
    IF (status /= 0) THEN
       CALL write_error('Error in nf90_put_var STD_STOKES')
       GOTO 999
    ENDIF

    IF (save_scats) THEN
       status = nf90_put_var(ncid,stokessingID,SUM_STOKES_SING(:,:,:,:),&
               & start=(/ iband,1,1,1,1 /),count=(/1,nmat,nx,ny,n_angles/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var SUM_STOKES_SING')
          GOTO 999
       ENDIF

       status = nf90_put_var(ncid,stokesmultID,SUM_STOKES_MULT(:,:,:,:),&
               & start=(/ iband,1,1,1,1 /),count=(/1,nmat,nx,ny,n_angles/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var SUM_STOKES_MULT')
          GOTO 999
       ENDIF

       status = nf90_put_var(ncid,contrib_stokesID,SUM_CONTRIB_STOKES_3D(:,:,:,0:nz_novac,:),&
               & start=(/ iband,1,1,1,1,1 /),count=(/1,nmat,nx,ny,nz_withsurf_novac,n_angles/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var SUM_CONTRIB_STOKES_3D')
          GOTO 999
       ENDIF

       status = nf90_put_var(ncid,singscat_contrib_stokesID,SUM_SINGSCAT_CONTRIB_STOKES_3D(:,:,:,0:nz_novac,:),&
               & start=(/ iband,1,1,1,1,1 /),count=(/1,nmat,nx,ny,nz_withsurf_novac,n_angles/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var SUM_SINGSCAT_CONTRIB_STOKES_3D')
          GOTO 999
       ENDIF       

       status = nf90_put_var(ncid,multscat_contrib_stokesID,SUM_MULTSCAT_CONTRIB_STOKES_3D(:,:,:,0:nz_novac,:),&
               & start=(/ iband,1,1,1,1,1 /),count=(/1,nmat,nx,ny,nz_withsurf_novac,n_angles/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var SUM_MULTSCAT_CONTRIB_STOKES_3D')
          GOTO 999
       ENDIF  
    ENDIF     

    IF (chatty) PRINT*,'Write DOP'
    status = nf90_put_var(ncid,dopID,DOP(:,:,:),start=(/ iband,1,1,1 /),count=(/1,nx,ny,n_angles/))
    IF (status /= 0) THEN
       PRINT*,'dopID=',dopID
       PRINT*,'ncid=',ncid
       PRINT*,'iband=',iband
       PRINT*,'DOP(1,1,1)=',DOP(1,1,1)
       CALL write_error('Error in nf90_put_var DOP')
       GOTO 999
    ENDIF

    IF (save_path) THEN
       IF (chatty) PRINT*,'Write paths'     
       status = nf90_put_var(ncid,amfID,AMF_layer(:,:,1:nz_novac,:),start=(/ iband,1,1,1,1 /),count=(/1,nx,ny,nz_novac,n_angles/))
       IF (status /= 0) THEN
          CALL write_error('Error in nf90_put_var AMF_layer')
          GOTO 999
       ENDIF
 
       !status = nf90_put_var(ncid,boxamfID,AMF_box(:,:,:,:,:,:),start=(/ iband,1,1,1,1,1,1 /),&
       !                      & count=(/1,nx,ny,nx,ny,nz_box,n_angles/))
       !IF (status /= 0) THEN
       !   CALL write_error('Error in nf90_put_var AMF_box')
       !   GOTO 999
       !ENDIF

       IF (chatty) PRINT*,'Writing paths finished'
    ENDIF


   !---------------------------------------------------------------------
   ! Increment number indicating we were in the first wav loop iteration
   !---------------------------------------------------------------------
   IF (chatty) PRINT*,'Increment number indicating we were in the first wav loop iteration'
   IF (outputnum.EQ.0) THEN
      outputnum = 1
   ENDIF

   !---------------------------------------------------------------------
   ! Close the output file    
   !---------------------------------------------------------------------
   IF (chatty) PRINT*,'Close the file'
   IF (outputnum.EQ.1.AND.iband.EQ.nbands) THEN
      Call my_nf90_close(ncid, status)
      IF (status /= 0) GOTO 999
   ENDIF

   !---------------------------------------------------------------------
   ! Print error statements if error occurred
   !---------------------------------------------------------------------
999 IF (status /= 0) THEN      
      CALL write_error(trim(nf90_strerror(status)))
      CALL write_error(err_str)
      CALL write_error('An error occured in Report_ncdf')      
   ENDIF

   RETURN

  END SUBROUTINE report_ncdf

END PROGRAM main
