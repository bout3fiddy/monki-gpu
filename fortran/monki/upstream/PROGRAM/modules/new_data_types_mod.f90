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

MODULE new_data_types
  
  IMPLICIT NONE
  
  TYPE new_scatterer_info_pol
     CHARACTER(len=80)             :: description
     INTEGER                       :: n_scat_angles ! Number of angles the phase function is tabulated at
     REAL*8,DIMENSION(:),POINTER     :: scat_angles   ! in deg
     REAL*8,DIMENSION(:),POINTER :: F11      ! F11 element of the phase function 
     Real*8,Dimension(:),Pointer :: F22
     Real*8,Dimension(:),Pointer :: F33
     Real*8,Dimension(:),Pointer :: F44
     Real*8,Dimension(:),Pointer :: F12
     Real*8,Dimension(:),Pointer :: F34
     Real*8,Dimension(:),Pointer :: CDF_bound
     Real*8,Pointer :: ssalb_aer
  END TYPE new_scatterer_info_pol

CONTAINS
  
!-----------------------------------------------------------------------------
! new_scatterer_info_pol
!-----------------------------------------------------------------------------
  SUBROUTINE allocate_and_initialize_new_scatterer_info_pol(in_str,n_scat_angles)

    TYPE(new_scatterer_info_pol),INTENT(INOUT) :: in_str
    INTEGER,INTENT(IN)                     :: n_scat_angles

    !Allocate
    ALLOCATE(in_str%scat_angles(1:n_scat_angles))
    ALLOCATE(in_str%CDF_bound(1:n_scat_angles))
    ALLOCATE(in_str%ssalb_aer)
    ALLOCATE(in_str%F11(1:n_scat_angles))
    ALLOCATE(in_str%F22(1:n_scat_angles))
    ALLOCATE(in_str%F33(1:n_scat_angles))
    ALLOCATE(in_str%F12(1:n_scat_angles))
    ALLOCATE(in_str%F34(1:n_scat_angles))
    ALLOCATE(in_str%F44(1:n_scat_angles))

    !Set all content equal to 0
    in_str%scat_angles = 0.D0
    in_str%CDF_bound(1:n_scat_angles) = 0.D0
    in_str%ssalb_aer = 0.D0
    in_str%F11(1:n_scat_angles) = 0.D0
    in_str%F22(1:n_scat_angles) = 0.D0
    in_str%F33(1:n_scat_angles) = 0.D0
    in_str%F12(1:n_scat_angles) = 0.D0
    in_str%F34(1:n_scat_angles) = 0.D0
    in_str%F44(1:n_scat_angles) = 0.D0
   
    in_str%n_scat_angles=n_scat_angles

  END SUBROUTINE allocate_and_initialize_new_scatterer_info_pol

END MODULE new_data_types
