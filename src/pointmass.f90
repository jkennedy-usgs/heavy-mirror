module POINTMASS
!***********************************************
! The pointmass module contains subroutines to
! calculate gravity from changes in water levels
! in model cells as if the mass change was
! concentrated at a point at the cell center.
!
! It is faster than the forsberg prism formula
! but the sensor must be far away from the model
! cell
!***********************************************
   implicit none
   doubleprecision, parameter, private :: gma = 6.673E-11
   doubleprecision, parameter, private :: rho = 1000.0
   doubleprecision, parameter, private :: ugal = 1e8

contains
   function PMGRAVITY(sy, hi, hf, posX, posY, posZ, xc, yc, &
                      delr, delc) result(grav)
      !*********************************************
      !Method to calculate relative gravity change
      !using the Pointmass approximation
      !
      !Parameters
      !   sy : real, (ncol, nrow)
      !       specific yield
      !   hi : real, (ncol, nrow)
      !       initial head
      !   hf : real, (ncol, nrow)
      !       head
      !   posX : real
      !       sensor x position
      !   posY : real
      !       sensor y position
      !   posZ : real
      !       sensor z position
      !   xc : real, (ncol, nrow)
      !       xcenters
      !   yc : real, (ncol, nrow)
      !       ycenters
      !   delr : real, (ncol, nrow)
      !       row cell size
      !   delc : real, (ncol, nrow)
      !       column cell size
      !
      !   Returns
      !       grav : real
      !       relative gravity change
      !*********************************************
      implicit none
      doubleprecision :: posX, posY, posZ
      doubleprecision :: x, y, z
      doubleprecision :: dh, sy, delr, delc
      doubleprecision :: xc, yc
      doubleprecision :: rf, grav, hi, hf
      x = posX - xc
      y = posY - yc
      z = posZ - hf
      dh = hf - hi
      rf = SQRT((x**2 + y**2 + z**2)**3)
      grav = ugal*gma*rho*sy*dh*delc*delr*(z/rf)

   end function
end module
