c	interpol.wsf				17-may-91
c
c	contains interpolation subroutines for linear, exp and
c	other interpolations
c
      subroutine linpol ( x, y, imax, xp, yp, jmax )
c     
c	Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp(1..jmax) (x-coordinates of the interpolated values)
c     
c     the interpolated y-values are returned in the array
c     
c     yp(1..jmax)
c     
      real  x(*), y(*), xp(*), yp(*)
      integer imax, jmax, i, j,imin
c     
      imin=1
c
c
      do 1 j = 1, jmax
         if (xp(j) .lt. x(1)) then 
            i=1
            yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
         else if (xp(j) .gt. x(imax)) then
            i=imax-1
            yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
         else if (xp(j) .eq. x(imax)) then
            yp(j)=y(imax)
         else
            do 2 i = imin, imax-1
               if (xp(j) .eq. x(i)) then
                  imin=i
                  yp(j)=y(i)
                  goto 1
               else if ((xp(j) .gt. x(i)).and.(xp(j) .le. x(i+1))) then
                  imin=i
c                  write(*,*) i,j,x(i+1),x(i),xp(j)
                  yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
                  goto 1
               endif
 2          continue
         endif
 1    continue
c     
      return
      end
c
      subroutine linpol_dp ( x, y, imax, xp, yp, jmax )
c     
c	Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp(1..jmax) (x-coordinates of the interpolated values)
c     
c     the interpolated y-values are returned in the array
c     
c     yp(1..jmax)
c     
      real*8  x(*), y(*), xp(*), yp(*)
      integer imax, jmax, i, j,imin
c     
      imin=1
c
c
      do 1 j = 1, jmax
         if (xp(j) .lt. x(1)) then 
            i=1
            yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
         else if (xp(j) .gt. x(imax)) then
            i=imax-1
            yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
         else if (xp(j) .eq. x(imax)) then
            yp(j)=y(imax)
         else
            do 2 i = imin, imax-1
               if (xp(j) .eq. x(i)) then
                  imin=i
                  yp(j)=y(i)
                  goto 1
               else if ((xp(j) .gt. x(i)).and.(xp(j) .le. x(i+1))) then
                  imin=i
                  !write(*,*) i,j,x(i+1),x(i),xp(j)
                  yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
                  goto 1
               endif
 2          continue
         endif
 1    continue
c     
      return
      end


      subroutine linpol_dp_single ( x, y, imax, xp, yp)
c     I
c       Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp (x-coordinate of the interpolated value)
c     
c     the interpolated y-value is returned in the number
c     
c     yp
c     
      real*8  x(*), y(*), xp, yp
      integer imax, i, imin
c     
      imin=1

      if (xp .lt. x(1)) then
         i=1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
      else if (xp .gt. x(imax)) then
         i=imax-1
         !yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i+1)) + y(i+1)
      else if (xp .eq. x(imax)) then
         yp=y(imax)
      else
         do 2 i = imin, imax-1
            if (xp .eq. x(i)) then
               !imin=i
               yp=y(i)
               goto 1
            !else if ((xp .gt. x(i)).and.(xp .le. x(i+1))) then
            else if ((xp .gt. x(i)).and.(xp .lt. x(i+1))) then
               !imin=i
               !write(*,*) i,j,x(i+1),x(i),xp(j)
               yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
               goto 1
            endif
 2       continue
      endif
c     
 1    return
      end


      subroutine linpol_dp_single_scat(x,y,imax,xp,yp)
c     I
c       Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp (x-coordinate of the interpolated value)
c     
c     the interpolated y-value is returned in the number
c     
c     yp
c     
      implicit none
      real*8  x(*), y(*), xp, yp
      integer imax, i, imin
      Real*8,parameter :: DegtoRad=1.7453292519943295D-002     ! Degrees to Radians 
c     
      imin=1

      if (xp .lt. x(1)) then
         i=1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
      else if (xp .gt. x(imax)) then
         i=imax-1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
      else if (xp .eq. x(imax)) then
         yp=y(imax)
      else

         !estimate where the index could be
         imin = INT(imax - DACOS(xp)/DegtoRad/180.D0*(imax))

         do 2 i = imin, imax-1
            if (xp .eq. x(i)) then
               imin=i
               yp=y(i)
               goto 1
            else if ((xp .gt. x(i)).and.(xp .le. x(i+1))) then
               imin=i
               !write(*,*) i,j,x(i+1),x(i),xp(j)
               yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
               goto 1
            endif
 2       continue            
         PRINT *,'ERROR in linpol_dp_single_scat:'
         PRINT *,'interpolation index not found!'
         PRINT *,'i=',i,'imin=',imin,'imax=',imax
         PRINT *,'xp=',xp
         READ(*,*)
      endif
c     
 1    return
      end


      subroutine linpol_dp_single1 ( x, y, imax, xp, yp, icase, imin)
c     I
c       Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp (x-coordinate of the interpolated value)
c     
c     the interpolated y-value is returned in the number
c     
c     yp
c     
      real*8  x(*), y(*), xp, yp
      integer imax, i, imin, icase
c     
      imin=1
      icase=0

c      do 1 j = 1, jmax
      if (xp .lt. x(1)) then
         i=1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
         icase=1
      else if (xp .gt. x(imax)) then
         i=imax-1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
         icase=2
      else if (xp .eq. x(imax)) then
         yp=y(imax)
         icase=3
      else
         do 2 i = imin, imax-1
            if (xp .eq. x(i)) then
               imin=i
               yp=y(i)
               icase=4
               goto 1
            else if ((xp .gt. x(i)).and.(xp .le. x(i+1))) then
               imin=i
               !write(*,*) i,j,x(i+1),x(i),xp(j)
               yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
               icase=5
               goto 1
            endif
 2       continue
      endif
c 1    continue
c     
 1    return
      end


      subroutine linpol_dp_single1_scat(x,y,imax,xp,yp,icase,imin)
c     I
c       Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp (x-coordinate of the interpolated value)
c     
c     the interpolated y-value is returned in the number
c     
c     yp
c     
      implicit none

      real*8 :: x(*), y(*), xp, yp
      integer :: imax, i, imin, icase
      Real*8,parameter :: DegtoRad=1.7453292519943295D-002     ! Degrees to Radians      
c     
      imin=1
      icase=0

c      do 1 j = 1, jmax
      if (xp .lt. x(1)) then
         i=1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
         icase=1
      else if (xp .gt. x(imax)) then
         i=imax-1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
         icase=2
      else if (xp .eq. x(imax)) then
         yp=y(imax)
         icase=3
      else

         !estimate where the index could be
         imin = MAX(1,INT(imax - DACOS(xp)/DegtoRad/180.D0*(imax)))

         do 2 i = imin, imax-1
            if (xp .eq. x(i)) then
               imin=i
               yp=y(i)
               icase=4
               goto 1
            else if ((xp .gt. x(i)).and.(xp .le. x(i+1))) then
               imin=i
               !write(*,*) i,j,x(i+1),x(i),xp(j)
               yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
               icase=5
               goto 1
            endif
 2       continue            
         PRINT *,'ERROR in linpol_dp_single1_scat:'
         PRINT *,'interpolation index not found!'
         PRINT *,'i=',i,'imin=',imin,'imax=',imax
         PRINT *,'xp=',xp               
         READ(*,*)
      endif
c 1    continue
c     
 1    return
      end


      subroutine linpol_dp_single2 ( x, y, imax, xp, yp, icase, imin)
c     I
c       Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp (x-coordinate of the interpolated value)
c     
c     the interpolated y-value is returned in the number
c     
c     yp
c     
      real*8  x(*), y(*), xp, yp
      integer imax, i, imin, icase
c     
      if (icase .eq. 1) then
         i=1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
      else if (icase .eq. 2) then
         i=imax-1
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
      else if (icase .eq. 3) then
         yp=y(imax)
      else if (icase .eq. 4) then
         i=imin
         yp=y(i)
      else if (icase .eq. 5) then
         i=imin
         yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
      else
         PRINT *,'ERROR in linpol_dp_single2: icase should be nonzero.'
         READ(*,*)
      endif
c     
      return
      end      



c
      subroutine linpol_sing ( x, y, imax, xp, yp)
c     
c	Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp(1) (x-coordinate of the interpolated values)
c     
c     the interpolated y-value are returned in the array
c     
c     yp(1)
c     
      real  x(*), y(*), xp, yp
      integer imax, jmax, i, j,imin
c     
      imin=1
c
c      write(*,*) x(1:imax)
c      write(*,*) '-------'
c      write(*,*) y(1:imax)
c
      jmax=1
c
      do 1 j = 1, jmax
         if (xp .lt. x(1)) then 
            i=1
            yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
         else if (xp .gt. x(imax)) then
            i=imax-1
            if (abs((x(i+1)-x(i))/(xp+x(i))).gt.1.0e-8) then
               yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
            else
               yp=y(i+1)
            endif
         else if (xp .eq. x(imax)) then
            yp=y(imax)
         else
            do 2 i = imin, imax-1
               if (xp .eq. x(i)) then
                  imin=i
                  yp=y(i)
                  goto 1
               else if ((xp .gt. x(i)).and.(xp .le. x(i+1))) then
                  imin=i
                  yp=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp-x(i)) + y(i)
                  goto 1
               endif
 2          continue
         endif
 1    continue
c     
      return
      end
c
c
c
      subroutine steppol ( x, y, imax, xp, yp, jmax )
c     
c	Does a step interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp(1..jmax) (x-coordinates of the interpolated values)
c     
c     the interpolated y-values are returned in the array
c     
c     yp(1..jmax)
c     
      real  x(*), y(*), xp(*), yp(*)
      integer imax, jmax, i, j,imin
c     
      imin=1
c
      do 1 j = 1, jmax
         if (xp(j) .le. x(1)) then 
            i=1
            yp(j)=y(i)
         else if (xp(j) .ge. x(imax)) then
            i=imax
            yp(j)=y(i)
         else
            do 2 i = imin+1, imax
               if ((xp(j) .le. x(i)).and.(xp(j) .gt. x(i-1))) then
                  imin=i-1
                  yp(j)=y(i)
                  goto 1
               endif
 2          continue
         endif
 1    continue
c     
      return
      end
c
c
c
c
      subroutine exppol ( x, y, imax, xp, yp, jmax )
c     
c     Does a  interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     A dependency y=exp(b*x) is assumed and so a linear interpolation
c     is done on yy=ln(y)=b*x. 
c     
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp(1..jmax) (x-coordinates of the interpolated values)
c     
c     the interpolated y-values are returned in the array
c     
c     yp(1..jmax)
c     
      real  x(*), y(*), xp(*), yp(*)
      integer imax, jmax, i, j,imin,idone
c     
      real  yip1,yi
c     
c     
      imin=1
      idone=-10
cd    write (*,*) 'imax=',imax,'jmax=',jmax
c     
      do 1 j = 1, jmax
c     
c     do extrapolation if xp(j) is smaller than x(1)
c     
         if (xp(j) .lt. x(1)) then 
            i=1
cd    write (*,*) 'xp(',j,')=',xp(j)
cd    write (*,*) 'extrapolating from x(',i,')=',x(i),
cd    +			    'and x(',i+1,')=',x(i+1)
cd    write (*,*) 'y(',i,')=',y(i),'and y(',i+1,')=',y(i+1)
            yip1=log(y(i+1))
            yi=log(y(i))
            yp(j)=exp((yip1-yi)/(x(i+1)-x(i))*(xp(j)-x(i)) + yi)
cd    write (*,*) 'gives yp(',j,')=',yp(j)
            goto 1
         endif
c     
c     do extrapolation if xp(j) is greater than x(imax)
c     
         if (xp(j) .gt. x(imax)) then
            i=imax-1
cd    write (*,*) 'xp(',j,')=',xp(j)
cd    write (*,*) 'extrapolating from x(',i,')=',x(i),
cd    +			    'and x(',i+1,')=',x(i+1)
cd    write (*,*) 'y(',i,')=',y(i),'and y(',i+1,')=',y(i+1)
            yip1=log(y(i+1))
            yi=log(y(i))
            yp(j)=exp((yip1-yi)/(x(i+1)-x(i))*(xp(j)-x(i)) + yi)
cd    write (*,*) 'gives yp(',j,')=',yp(j)
            goto 1
         endif
c     
c     check if xp(j) is x(imax)
c     
         if (xp(j) .eq. x(imax)) then
            yp(j)=y(imax)
            goto 1
         endif
c     
c     find out between which x(i) values the actual xp(i) lies
c     
         do 2 i = imin, imax-1
            if (xp(j) .eq. x(i)) then
	       imin=i
	       yp(j)=y(i)
	       goto 1
            else if ((xp(j) .gt. x(i)).and.(xp(j) .lt. x(i+1))) then
c     
c     check if we already calculated the logarithms
c     
	       if (i .eq. idone) goto 10
c     
c     check if we already calculated one of the logs
c     
	       if (i .eq. idone+1) then
                  yi=yip1
                  goto 5
	       endif
c     
c     otherwise calculate both logs
c     
	       yi=log(y(i))
 5             yip1=log(y(i+1))
c     
c     do the interpolation
c     
 10            imin=i
	       idone=i
cd    write (*,*) 'xp(',j,')=',xp(j)
cd    write (*,*) 'interpolating from x(',i,')=',x(i),
cd    +			    'and x(',i+1,')=',x(i+1)
cd    write (*,*) 'y(',i,')=',y(i),'and y(',i+1,')=',y(i+1)
	       yp(j)=exp((yip1-yi)/(x(i+1)-x(i))*(xp(j)-x(i)) + yi)
cd    write (*,*) 'gives yp(',j,')=',yp(j)
	       goto 1
            endif
 2       continue
 1    continue
c     
      return
      end



c
      subroutine exppol_sing( x, y, imax, xp, yp)
c
      implicit none
c     
c     Does a  interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     A dependency y=exp(b*x) is assumed and so a linear interpolation
c     is done on yy=ln(y)=b*x. 
c     
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp (x-coordinate of the interpolated values)
c     
c     the interpolated y-values are returned in the value
c     
c     yp
c     
      real  x(*), y(*)
      Real  xp, yp
      integer imax, i,imin,idone
c     
      real  yip1,yi
c
c     
      yip1=y(1)
c     
      imin=1
      idone=-100
c
c     
c     do extrapolation if xp(j) is smaller than x(1)
c     
      if (xp .lt. x(1)) then 
         i=1
         yip1=log(y(i+1))
         yi=log(y(i))
         yp=exp((yip1-yi)/(x(i+1)-x(i))*(xp-x(i)) + yi)
         goto 1
      endif
c     
c     do extrapolation if xp(j) is greater than x(imax)
c     
      if (xp .gt. x(imax)) then
         i=imax-1
         yip1=log(y(i+1))
         yi=log(y(i))
         yp=exp((yip1-yi)/(x(i+1)-x(i))*(xp-x(i)) + yi)
         goto 1
      endif
c     
c     check if xp(j) is x(imax)
c     
      if (xp .eq. x(imax)) then
         yp=y(imax)
         goto 1
      endif
c     
c     find out between which x(i) values the actual xp(i) lies
c     
      imin=1
      do 2 i = imin, imax-1
         if (xp .eq. x(i)) then
            imin=i
            yp=y(i)
            goto 1
         else if ((xp .gt. x(i)).and.(xp .lt. x(i+1))) then
c     
c     check if we already calculated the logarithms
c     
            if (i .eq. idone) goto 10
c     
c     check if we already calculated one of the logs
c     
            if (i .eq. idone+1) then
               yi=yip1
               goto 5
            endif
c     
c     otherwise calculate both logs
c     
            yi=log(y(i))
 5          yip1=log(y(i+1))
c     
c     do the interpolation
c     
 10         imin=i
            idone=i
            yp=exp((yip1-yi)/(x(i+1)-x(i))*(xp-x(i)) + yi)
            goto 1
         endif
 2    continue
 1    return
      end










      subroutine linpol_nord ( x, y, imax, xp, yp, jmax )
c     
c	Does a linear interpolation on the data given by the arrays
c     
c     x(1..imax) (x-coordinate) and
c     y(1..imax) (y-coordinate)
c     
c     Data are interpolated at the x-values given by the array
c     
c     xp(1..jmax) (x-coordinates of the interpolated values)
c     
c     the interpolated y-values are returned in the array
c     
c----------------------------------------------------
c     xp does not have to be in ascending order !!!!
c----------------------------------------------------
c
c     yp(1..jmax)
c     
      real  x(*), y(*), xp(*), yp(*)
      integer imax, jmax, i, j,imin
c     
      imin=1
c
      do 1 j = 1, jmax
         if (xp(j) .lt. x(1)) then 
            i=1
            yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
            goto 1
         endif
         if (xp(j) .gt. x(imax)) then
            i=imax-1
            yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
            goto 1
         endif
         if (xp(j) .eq. x(imax)) then
            yp(j)=y(imax)
            goto 1
         endif
         do 2 i = 1, imax-1
            if (xp(j) .eq. x(i)) then
	       yp(j)=y(i)
	       goto 1
            else if ((xp(j) .gt. x(i)).and.(xp(j) .le. x(i+1))) then
	       yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
	       goto 1
            endif
 2       continue
 1    continue
c     
      return
      end
















