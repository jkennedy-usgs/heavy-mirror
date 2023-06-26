module FORSBERG
!***********************************************
! The forsberg module contains subroutines to
! calculate gravity from changes in water levels
! supplied x, y, z information
!***********************************************
   implicit none

   doubleprecision, parameter, private :: gma = 6.673d-11
   doubleprecision, parameter, private :: rho = 1000.0d0
   doubleprecision, parameter, private :: ugal = 1d8

contains
   function GRAVITY(sy, hi, hf, posX, posY, posZ, xc, yc, delr, delc) result(grav)
      !************************************************
      !Method to calculate relative gravity change
      !
      !Parameters
      !   sy : real
      !       specific yield
      !   hi : real
      !       initial head
      !   hf : real
      !       final head
      !   posX : real
      !       x-position of gravity measurement
      !   posY : real
      !       y-position of gravity measurement
      !   posZ : real
      !       z-position of gravity measurement
      !   xc : real
      !       x-position of cell center
      !   yc : real
      !       y-position of cell center
      !   delr : real
      !       cells delr from MODFLOW
      !   delc : real
      !       cells delc from MODFLOW
      !
      !Returns
      !   grav : real
      !       relative gravity change
      !*************************************************
      implicit none
      doubleprecision :: sy, posX, posY, posZ, delr, delc
      doubleprecision :: dx0, dx1, dy0, dy1, dz0, dz1, xc, yc
      doubleprecision :: grav, hf, hi
      doubleprecision, dimension(8) :: x, y, z
      doubleprecision, dimension(8) :: rf
      doubleprecision, dimension(8) :: eigen
      doubleprecision, dimension(8) :: gvarr
      doubleprecision :: temp_sum

      dx0 = (xc - (delc/2.d0)) - posX
      dx1 = (xc + (delc/2.d0)) - posX
      dy0 = (yc - (delr/2.d0)) - posY
      dy1 = (yc + (delr/2.d0)) - posY
      dz0 = hf - posZ
      dz1 = hi - posZ

      x = [dx0, dx0, dx0, dx0, dx1, dx1, dx1, dx1]
      y = [dy0, dy0, dy1, dy1, dy0, dy0, dy1, dy1]
      z = [dz0, dz1, dz0, dz1, dz0, dz1, dz0, dz1]
      eigen = [1.d0, -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0]

      rf = SQRT(x**2.d0 + y**2.d0 + z**2.d0)
      grav = 0
      IF (all(rf .ne. 0)) THEN
         IF (all(z .ne. 0)) THEN
            gvarr = (eigen*(x*LOG(y + rf) + y* &
                            LOG(x + rf) - z*ATAN(x*y/z/rf)))
            temp_sum = SUM(gvarr)
            grav = gma*sy*rho*ugal*SUM(gvarr)
         END IF
      END IF
      ! open(61,file='debug.txt',action='write',position='append')
      ! write(61,  *) x, y, z, temp_sum

   end function

end module
