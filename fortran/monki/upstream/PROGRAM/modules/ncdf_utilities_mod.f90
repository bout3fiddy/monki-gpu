  !! *MODULE* ncdf_utilities
  !!
  !! @version ecsim 1.2.4
  !!
  !! *SRC_FILE*
  !!
  !! models/common/src/ncdf_utilities_mod.f90
  !!
  !! *LAST CHANGES*
  !!
  !!  -Oct 08, 2012: D.D. Added 4 and 5 d put int routines
  !!  -Oct 01, 2012: G.J. Added 99 as input possibility to status in open_ncdf. This
  !!                      will ensure that no error message is produced
  !!  -Sep 24, 2012: G.J. Added definition of status (integer) to all subroutines
  !!  -May 24, 2011: D.D. Changed errors to warnings in some of the routines
  !!  -Mar 10, 2008: D.D. Modified the header and Copyright info
  !!
  !! *DESCRIPTION* 
  !!
  !! Module that contains a number of subroutines to support the handling of NetCDF data. 
  !!
  ! 
  !  Copyright KNMI for ESA
  !
  !  This source code is free software: you can redistribute it and/or modify
  !  it under the terms of the GNU Lesser General Public License as published by
  !  the Free Software Foundation, either version 3 of the License, or any later version.
  !
  !  This program is distributed in the hope that it will be useful,
  !  but WITHOUT ANY WARRANTY; without even the implied warranty of
  !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  !  GNU Lesser General Public License for more details.
  !
  !  You should have received a copy of the GNU Lesser General Public License
  !  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  !
Module ncdf_utilities
! ---------------------
! Version March 2007
! Author: Robert Voors
! ---------------------
Use write_messages
Use typeSizes
Use netcdf
!
Interface my_nf90_put_att
!   module procedure my_nf90_put_att_text, my_nf90_put_att_text_1d, &
   Module Procedure my_nf90_put_att_text,  &       ! No 1-D text...
                    my_nf90_put_att_int,  my_nf90_put_att_int_1d,  &
                    my_nf90_put_att_real, my_nf90_put_att_real_1d, &
                    my_nf90_put_att_dbl,  my_nf90_put_att_dbl_1d
End Interface ! my_nf90_put_att
!
Interface my_nf90_put_var
  Module Procedure my_nf90_put_var_int_0d,  my_nf90_put_var_int_1d,  my_nf90_put_var_int_2d,  &
       & my_nf90_put_var_int_4d,  my_nf90_put_var_int_5d,  &
       &my_nf90_put_var_int_3d,  my_nf90_put_var_real_0d, my_nf90_put_var_real_1d, &
       &my_nf90_put_var_real_2d, my_nf90_put_var_real_3d, my_nf90_put_var_dbl_0d,  &
       &my_nf90_put_var_dbl_1d,  my_nf90_put_var_dbl_2d,  my_nf90_put_var_dbl_3d
End Interface ! my_nf90_put_var
!
Interface my_nf90_get_var
  Module Procedure my_nf90_get_var_int_0d,  my_nf90_get_var_int_1d,  my_nf90_get_var_int_2d,  &
                  &my_nf90_get_var_int_3d,  my_nf90_get_var_real_0d, my_nf90_get_var_real_1d, &
                  &my_nf90_get_var_real_2d, my_nf90_get_var_real_3d, my_nf90_get_var_dbl_0d,  &
                  &my_nf90_get_var_dbl_1d,  my_nf90_get_var_dbl_2d,  my_nf90_get_var_dbl_3d, &
                  &my_nf90_get_var_real_4d
End Interface ! my_nf90_get_var
!
!
Contains
!
Subroutine my_nf90_get_var_int_0d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Intent(out)                        :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_int_0d
!
!
Subroutine my_nf90_get_var_int_1d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Dimension(:), Intent(out)          :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_int_1d
!
!
Subroutine my_nf90_get_var_int_2d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Dimension(:,:), Intent(out)        :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_int_2d
!
!
Subroutine my_nf90_get_var_int_3d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Dimension(:,:,:), Intent(out)      :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
   Integer                                      :: status
 !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_int_3d
!
!
Subroutine my_nf90_get_var_real_0d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Intent(out)                           :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_real_0d
!
!
Subroutine my_nf90_get_var_real_1d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Dimension(:), Intent(out)             :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_real_1d
!
!
Subroutine my_nf90_get_var_real_2d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Dimension(:,:), Intent(out)           :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_real_2d
!
!
Subroutine my_nf90_get_var_real_3d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Dimension(:,:,:), Intent(out)         :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_real_3d
!
!
Subroutine my_nf90_get_var_real_4d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Dimension(:,:,:,:), Intent(out)       :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_real_4d
!
!
Subroutine my_nf90_get_var_dbl_0d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real*8, Intent(out)                         :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_dbl_0d
!
!
Subroutine my_nf90_get_var_dbl_1d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real*8, Dimension(:), Intent(out)           :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_dbl_1d
!
!
Subroutine my_nf90_get_var_dbl_2d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real*8, Dimension(:,:), Intent(out)         :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_dbl_2d
!
!
Subroutine my_nf90_get_var_dbl_3d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real*8, Dimension(:,:,:), Intent(out)       :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_get_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = ' in nf90_get_var: '
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_get_var_dbl_3d
!






! **************************************************8
! ***************************************************


!
Subroutine my_nf90_put_var_int_0d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Intent(in)                         :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
!  integer, dimension(:), optional :: start, count, stride, map
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_int_0d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_int_0d
!
!
Subroutine my_nf90_put_var_int_1d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Dimension(:), Intent(in)           :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_int_1d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_int_1d
!
!
Subroutine my_nf90_put_var_int_2d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Dimension(:,:), Intent(in)         :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_int_2d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_int_2d
!
!
Subroutine my_nf90_put_var_int_3d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Dimension(:,:,:), Intent(in)       :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_int_3d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_int_3d
!
Subroutine my_nf90_put_var_int_4d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Dimension(:,:,:,:), Intent(in)       :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_int_4d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_int_4d
!
Subroutine my_nf90_put_var_int_5d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Integer, Dimension(:,:,:,:,:), Intent(in)       :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_int_5d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_int_5d
!
!
Subroutine my_nf90_put_var_real_0d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Intent(in)                            :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_real_0d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_real_0d
!
!
Subroutine my_nf90_put_var_real_1d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Dimension(:), Intent(in)              :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_real_1d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_real_1d
!
!
Subroutine my_nf90_put_var_real_2d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Dimension(:,:), Intent(in)            :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_real_2d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_real_2d
!
!
Subroutine my_nf90_put_var_real_3d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real, Dimension(:,:,:), Intent(in)          :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_real_3d:'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_real_3d
!

!
Subroutine my_nf90_put_var_dbl_0d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real*8, Intent(in)                          :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_dbl_0d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_dbl_0d
!
!
Subroutine my_nf90_put_var_dbl_1d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real*8, Dimension(:), Intent(in)            :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_dbl_1d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_dbl_1d
!
!
Subroutine my_nf90_put_var_dbl_2d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real*8, Dimension(:,:), Intent(in)          :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_dbl_2d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_dbl_2d
!
!
Subroutine my_nf90_put_var_dbl_3d(ncid, varid, values, return_status) 
  Integer, Intent(in)                         :: ncid, varid
  Real*8, Dimension(:,:,:), Intent(in)        :: values
  Integer, Intent(out)                        :: return_status
  Character(len=nf90_max_name)                :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_var(ncid, varid, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_put_var_dbl_3d: '
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_var_dbl_3d
!

! ********************************************
! ********************************************

!
Subroutine my_nf90_def_var(ncid, name, Type, dimIDs, varID, return_status)
  Integer, Intent(in)               :: ncid, Type
  Character(len=*), Intent(in)      :: name
  Integer, Dimension(:), Intent(in) :: dimIDs 
  Integer, Intent(out)              :: return_status, varID
  Character(len=nf90_max_name)      :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_def_var(ncid, name, Type, dimIDs, varID)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_def_var: '//name
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_def_var
!
!
Subroutine my_nf90_def_dim(ncid, name, dimsize, dimID, return_status)
  Integer, Intent(in)          :: ncid, dimsize
  Character(len=*), Intent(in) :: name
  Integer, Intent(out)         :: return_status, dimID
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_def_dim(ncid, Trim(name), dimsize, dimID)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_def_dim: '//name
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_def_dim
!
!
Subroutine my_nf90_create(name, cmode, ncid, return_status)
  Integer, Intent(in)          :: cmode
  Character(len=*), Intent(in) :: name
  Integer, Intent(out)         :: ncid, return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_create(name, cmode, ncid)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_create: '//name
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_create
!
!
Subroutine my_nf90_enddef(ncid, return_status)
  Integer, Intent(in)          :: ncid
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_enddef(ncid)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_enddef'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_enddef
!
!
Subroutine my_nf90_close(ncid, return_status)
  Integer, Intent(in)          :: ncid
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_close(ncid)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_enddef'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_close
!
!
Subroutine my_nf90_open(name, mode, ncid, return_status)
  Character(len=*), Intent(in) :: name
  Integer, Intent(in)          :: mode
  Integer, Intent(out)         :: ncid, return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                      :: status
  !
  if (return_status.ne.99) then
     return_status = 0
  endif
  status = nf90_open(name, mode, ncid)

  If (status /= 0) Then 
     err_msg  = 'Error in nf90_open'
     if (return_status.eq.0) Call write_error(err_msg)
     return_status = 1
  ELSE
     return_status = 0
  Endif
  
  Return
End Subroutine my_nf90_open
!
!
Subroutine my_nf90_inquire(ncid, return_status, nDims, nVars, nGlAtts, unlimitedId, formatNum)
  Integer, Intent(in)            :: ncid
  Integer, Optional, Intent(out) :: nDims, nVars, nGlAtts, unlimitedId, formatNum
  Integer, Intent(out)           :: return_status
  Character(len=nf90_max_name)   :: err_msg
   Integer                                      :: status
 !
  return_status = 0
  status = nf90_inquire(ncid, nDims, nVars, nGlAtts, unlimitedId, formatNum)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_inquire'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_inquire
!
!
Subroutine my_nf90_inquire_dimension(ncid, dimId, return_status, name, len)
  Integer, Intent(in)                     :: ncid, dimid
  Character(len=*), Optional, Intent(out) :: name
  Integer, Optional, Intent(out)          :: len
  Integer, Intent(out)                    :: return_status
  Character(len=nf90_max_name)            :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_inquire_dimension(ncid, dimId, name, len)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_inquire_dimension'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_inquire_dimension
!
!
Subroutine my_nf90_inquire_variable(ncid, varId, return_status, name, xtype, ndims, dimids, nAtts)
  Integer, Intent(in)                          :: ncid, varid
  Character(len=*), Optional, Intent(out)      :: name
  Integer, Optional, Intent(out)               :: xtype, ndims, nAtts
  Integer, Dimension(:), Optional, Intent(out) :: dimids 
  Integer, Intent(out)                         :: return_status
  Character(len=nf90_max_name)                 :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_inquire_variable(ncid, varId, name, xtype, ndims, dimids, nAtts) 
  If (status /= 0) Then 
     err_msg  = 'in nf90_inquire_variable: '//name
     Call write_warning(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_inquire_variable
!
!
Subroutine my_nf90_inq_varid(ncid, return_status, name, varid)
  Integer, Intent(in)                          :: ncid
  Character(len=*), Optional, Intent(in)       :: name
  Integer, Intent(out)                         :: return_status
  Integer, Optional, Intent(out)               :: varid
  Character(len=nf90_max_name)                 :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_inq_varid(ncid, name, varid)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_inq_varid: '//name
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_inq_varid
!

!
! PUT ATTRIBUTES
!
Subroutine my_nf90_put_att_text(ncid, varid, name, values, return_status)
  Integer, Intent(in)          :: ncid, varid
  Character(len=*), Intent(in) :: name
  Character(len=*), Intent(in) :: values
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  !print*, ncid, varid, name, values
  return_status = 0
  status = nf90_put_att(ncid, varid, name, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_putt_att_text'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_att_text
!
!
Subroutine my_nf90_put_att_int(ncid, varid, name, values, return_status)
  Integer, Intent(in)          :: ncid, varid
  Character(len=*), Intent(in) :: name
  Integer, Intent(in)          :: values
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_att(ncid, varid, name, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_putt_att_int'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_att_int
!
!
Subroutine my_nf90_put_att_real(ncid, varid, name, values, return_status)
  Integer, Intent(in)          :: ncid, varid
  Character(len=*), Intent(in) :: name
  Real, Intent(in)             :: values
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_att(ncid, varid, name, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_putt_att_real'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_att_real
!
!
Subroutine my_nf90_put_att_dbl(ncid, varid, name, values, return_status)
  Integer, Intent(in)          :: ncid, varid
  Character(len=*), Intent(in) :: name
  Real*8, Intent(in)           :: values
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_att(ncid, varid, name, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_putt_att_dbl'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_att_dbl
!
! 1-D
!
!Subroutine my_nf90_put_att_text_1d(ncid, varid, name, values, return_status)
!  integer, intent(in)          :: ncid, varid
!  character(len=*), intent(in) :: name
!  character(len=*), dimension(:), intent(in) :: values
!  integer, intent(out)         :: return_status
!  character(len=nf90_max_name) :: err_msg
!  !
!  return_status = 0
!  status = nf90_put_att(ncid, varid, name, values)
!  if (status /= 0) then 
!     err_msg  = 'Error in nf90_putt_att_text_1d'
!     call write_error(err_msg)
!     return_status = 1
!  endif
!  return
!End Subroutine my_nf90_put_att_text_1d
!
!
Subroutine my_nf90_put_att_int_1d(ncid, varid, name, values, return_status)
  Integer, Intent(in)          :: ncid, varid
  Character(len=*), Intent(in) :: name
  Integer, Dimension(:), Intent(in)          :: values
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_att(ncid, varid, name, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_putt_att_int_1d'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_att_int_1d
!
!
Subroutine my_nf90_put_att_real_1d(ncid, varid, name, values, return_status)
  Integer, Intent(in)          :: ncid, varid
  Character(len=*), Intent(in) :: name
  Real, Dimension(:), Intent(in)             :: values
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_att(ncid, varid, name, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_putt_att_real_1d'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_att_real_1d
!
!
Subroutine my_nf90_put_att_dbl_1d(ncid, varid, name, values, return_status)
  Integer, Intent(in)          :: ncid, varid
  Character(len=*), Intent(in) :: name
  Real*8, Dimension(:), Intent(in)           :: values
  Integer, Intent(out)         :: return_status
  Character(len=nf90_max_name) :: err_msg
  Integer                                      :: status
  !
  return_status = 0
  status = nf90_put_att(ncid, varid, name, values)
  If (status /= 0) Then 
     err_msg  = 'Error in nf90_putt_att_dbl_1d'
     Call write_error(err_msg)
     return_status = 1
  Endif
  Return
End Subroutine my_nf90_put_att_dbl_1d
!

!
Subroutine add_attributes(ncid,varID,long_name,units,plot_title, return_status)
  !
  Integer,Intent(in)          :: ncid,varID
  Character(len=*),Intent(in) :: long_name
  Character(len=*),Intent(in) :: units
  Character(len=*),Intent(in) :: plot_title
  integer, intent(out)        :: return_status
  !
  Integer :: status
  Character(len=200)          :: err_msg
  !
  status=nf90_put_att(ncid,varID,"long_name",long_name)
  If (status /= 0) Then
     err_msg = 'Error in nf90_put_att: '//Trim(Adjustl(long_name))
     Goto 999
  Endif
  !
  status=nf90_put_att(ncid,varID,"units",units)
  If (status /= 0) Then
     err_msg = 'Error in nf90_put_att: '//Trim(Adjustl(units))
     Goto 999
  Endif
  !
  status=nf90_put_att(ncid,varID,"plot_title",plot_title)
  If (status /= 0) Then
     err_msg = 'Error in nf90_put_att: '//Trim(Adjustl(plot_title))
     Goto 999
  Endif
  !
999 If (status.Ne.0) Then 
     Call write_error(err_msg)
     return_status = status
  Endif
  !
End Subroutine add_attributes
!
End Module ncdf_utilities
