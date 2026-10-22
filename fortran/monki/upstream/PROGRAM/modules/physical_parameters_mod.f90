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

module physical_parameters
  
  real*8,public,parameter :: pi=3.141592653589793D0
  Real*8,public,parameter :: DegtoRad=1.7453292519943295D-002 ! Degrees to Radians 
  Real*8,public,parameter :: avogad=6.022137D23    !DAK KNMI
  !Real*8,public,parameter :: loschmidt=2.54743D25 !DAP
  Real*8,public,parameter :: loschmidt=2.5468993D25 !DAK KNMI, number of air molecules per m3 at 288.15K and 1013.25 hPa
  REAL*8,public,parameter :: gasconst=8.314510D0 !Universal gas constant (in N m deg^-1 mol^-1):  (IUPAP, 1987)

end module physical_parameters



